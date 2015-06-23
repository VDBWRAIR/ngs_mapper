from bwa.bwa import BWAMem, index_ref, which_bwa, compile_refs
from ngs_mapper.data import reads_by_plat
from ngs_mapper.reads import compile_reads
import ngs_mapper.bam

import os
import sys
from os.path import *
import log
import tempfile
import shutil
from subprocess import PIPE

logger = log.setup_logger(__name__, log.get_config())

# For bwa errors
class BWAError(Exception): pass

def main():
    '''
        Compiles and runs everything

        @returns path to the final bam file which will be dictated by the --ouput arg value
    '''
    args = parse_args( sys.argv[1:] )
    tdir = join(dirname(args.output), 'bwa')
    
    # Compile together all the reads into a list
    preads = reads_by_plat( args.reads )
    logger.debug( "Reads parsed by platform: {0}".format(preads) )
    reads = []
    for plat in args.platforms:
        if plat in preads:
            reads += preads[plat]
    # Creates reads/F.fq, reads/R.fq, reads/NP.fq
    readdir = join(tdir,'reads')
    os.makedirs( readdir )
    reads = compile_reads( reads, readdir )
    if not reads:
        raise Exception( "Somehow no reads were compiled" )

    if os.path.isdir( args.reference ):
        cwd = os.getcwd()
        os.chdir( tdir )
        ref = join( tdir, compile_refs( args.reference ) )
        os.chdir( cwd )
    else:
        ref = args.reference

    # Keeps track so we know to merge bams later if it is 3
    merge = 0

    if reads['F'] is not None:
        merge += 1
        pairedsai = bwa_mem( reads['F'], reads['R'], ref, join(tdir, 'paired.sai'), t=args.threads )
        if isinstance(pairedsai,int):
            raise BWAError("There was an error running bwa")
        pairedbam = ngs_mapper.bam.sortbam( ngs_mapper.bam.samtobam( pairedsai, PIPE ), join(tdir, 'paired.bam') )
        #bam.indexbam( pairedbam )

    if reads['NP'] is not None:
        merge += 2
        nonpairedsai = bwa_mem( reads['NP'], ref=ref, output=join(tdir, 'nonpaired.sai'), t=args.threads )
        if isinstance(nonpairedsai,int):
            raise BWAError("There was an error running bwa")
        nonpairedbam = ngs_mapper.bam.sortbam( ngs_mapper.bam.samtobam( nonpairedsai, PIPE ), join(tdir, 'nonpaired.bam') )
        #bam.indexbam( nonpairedbam )

    # Now decide if any merging needs to happen
    bampath = args.output
    if merge == 3:
        ngs_mapper.bam.mergebams( [pairedbam, nonpairedbam], args.output )
    elif merge == 1:
        logger.debug( "Paired only. Moving result file {0} to {1}".format(pairedbam, bampath) )
        shutil.move( pairedbam, bampath )
    elif merge == 2:
        logger.debug( "Paired only. Moving result file {0} to {1}".format(nonpairedbam, bampath) )
        shutil.move( nonpairedbam, bampath )
    else:
        raise Exception( "Somehow no reads were compiled" )

    # Index the resulting bam
    ngs_mapper.bam.indexbam( bampath )

    if not args.keep_temp:
        shutil.rmtree( tdir )
    else:
        logger.info( "Keeping temporary directory {0}. You will probably want to delete it yourself or move it".format(tdir) )

def parse_args( args=sys.argv[1:] ):
    '''
        Uses argparse to parse arguments

        @params args - List of command line arguments to parse

        @returns Namespace object with parsed args in it(Same as ArgumentParser.parse_args)
    '''
    import argparse

    from ngs_mapper import config
    conf_parser, args, config, configfile = config.get_config_argparse(args)
    defaults = config['run_bwa_on_samplename']

    parser = argparse.ArgumentParser(
        description='Runs the bwa mem argument on a given set of reads and references for the given platform\'s reads',
        epilog='You should consider this script an autonomous bwa operation. That is, it will select the reads for ' \
            'the platforms you select and it will pull those reads onto the current host. It will ' \
            'compile references if they are multiple ones in a directory you select onto the local computer. Then it ' \
            'will map any mated reads against to those refs and also map the nonpaired reads against that ref in a separate ' \
            'call. Then when both are finished it will convert to bam/sort/index/merge/reindex the results. It attempts all of ' \
            'this inside of the /dev/shm filesystem which should be very fast. If /dev/shm cannot be used then /tmp will be used. '\
            'If you want the temporary files that are created to stay then you can use the --keep-temp argument',
        parents=[conf_parser]
    )
    
    parser.add_argument(
        dest='reads',
        help='Location to a directory containing fastq read files'
    )

    parser.add_argument(
        dest='reference',
        help='Location of the reference to map to. Should be a fasta file or directory of fasta files that will be compiled together.'
    )

    parser.add_argument(
        '--platforms',
        dest='platforms',
        nargs='+',
        choices=defaults['platforms']['choices'],
        default=defaults['platforms']['default'],
        help=defaults['platforms']['help']
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        default=defaults['output']['default'],
        help=defaults['output']['default']
    )

    parser.add_argument(
        '--keep-temp',
        dest='keep_temp',
        action='store_true',
        default=defaults['keep_temp']['default'],
        help=defaults['keep_temp']['help']
    )

    parser.add_argument(
        '-t',
        dest='threads',
        default=defaults['threads']['default'],
        type=int,
        help=defaults['threads']['help']
    )

    return parser.parse_args( args )

class InvalidReference(Exception): pass

def bwa_mem( read1, mate=None, ref=None, output='bwa.sai', **kwargs ):
    '''
        Runs the bwa mem algorithm on read1 against ref. If mate is given then run that file with the read1 file
        so paired alignment is done.

        TODO:
            bwa_path should be an option to specify where the executable is

        @param read1 - File path to read
        @param mate - Mate file path
        @param ref - Reference file path or directory of references
        @param output - The output destination

        @returns the output path if sucessful or -1 if something went wrong
    '''
    if os.path.isdir( ref ):
        # Compile ref directory
        logger.debug( "Compiling references inside of {0}".format(ref) )
        ref = compile_refs( ref )
        logger.info( "Refs are all compiled into {0}".format(ref) )

    # First, make sure the reference is indexed
    logger.debug( "Ensuring {0} is indexed".format(ref) )
    if not index_ref(ref):
        raise InvalidReference("{0} cannot be indexed by bwa")

    # Setup BWA Mem
    mem = None
    if mate:
        mem = BWAMem( ref, read1, mate, bwa_path=which_bwa(), **kwargs )
    else:
        mem = BWAMem( ref, read1, bwa_path=which_bwa(), **kwargs )

    ret = mem.run( output )
    print "Ret: " + str(ret)
    if ret != 0:
        return ret
    else:
        return output
