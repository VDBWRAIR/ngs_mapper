from bwa.bwa import BWAMem, index_ref, which_bwa, compile_refs
from data import reads_by_plat
from reads import compile_reads
import bam

import os
import sys
from os.path import *
import logging
import tempfile
import shutil
from subprocess import PIPE

log = logging.getLogger(__name__)

def main():
    '''
        Compiles and runs everything

        @returns path to the final bam file which will be dictated by the --ouput arg value
    '''
    args = parse_args( sys.argv[1:] )
    # Lets try shared memory
    # This is an untested code path right now
    # Probably scary stuff happens if /dev/shm fills up, but whatev
    if os.path.exists( '/dev/shm' ):
        tdir = tempfile.mkdtemp(prefix='mapbwa',dir='/dev/shm')
    else:
        tdir = tempfile.mkdtemp(prefix='mapbwa')
    
    # Compile together all the reads into a list
    preads = reads_by_plat( args.reads )
    reads = []
    for plat in args.platforms:
        if plat in preads:
            reads += preads[plat]
    # Creates reads/F.fq, reads/R.fq, reads/NP.fq
    readdir = join(tdir,'reads')
    os.mkdir( readdir )
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
        pairedsai = bwa_mem( reads['F'], reads['R'], ref, join(tdir, 'paired.sai') )
        pairedbam = bam.sortbam( bam.samtobam( pairedsai, PIPE ), join(tdir, 'paired.bam') )
        bam.indexbam( pairedbam )

    if reads['NP'] is not None:
        merge += 2
        nonpairedsai = bwa_mem( reads['NP'], ref=ref, output=join(tdir, 'nonpaired.sai') )
        nonpairedbam = bam.sortbam( bam.samtobam( nonpairedsai, PIPE ), join(tdir, 'nonpaired.bam') )
        bam.indexbam( nonpairedbam )

    # Now decide if any merging needs to happen
    bampath = args.output
    if merge == 3:
        bam.mergebams( [pairedbam, nonpairedbam], args.output )
        bam.indexbam( bampath )
    elif merge == 1:
        log.debug( "Paired only. Moving result file {} to {}".format(pairedbam, bampath) )
        shutil.move( pairedbam, bampath )
    elif merge == 2:
        log.debug( "Paired only. Moving result file {} to {}".format(nonpairedbam, bampath) )
        shutil.move( nonpairedbam, bampath )
    else:
        raise Exception( "Somehow no reads were compiled" )

    if not args.keep_temp:
        shutil.rmtree( tdir )
    else:
        log.info( "Keeping temporary directory {}. You will probably want to delete it yourself or move it".format(tdir) )

    return bampath

def parse_args( args=sys.argv[1:] ):
    '''
        Uses argparse to parse arguments

        @params args - List of command line arguments to parse

        @returns Namespace object with parsed args in it(Same as ArgumentParser.parse_args)
    '''
    import argparse
    parser = argparse.ArgumentParser(
        description='Runs the bwa mem argument on a given set of reads and references for the given platform\'s reads',
        epilog='You should consider this script an autonomous bwa operation. That is, it will select the reads for ' \
            'the platforms you select(MiSeq,Sanger only) and it will pull those reads onto the current host. It will ' \
            'compile references if they are multiple ones in a directory you select onto the local computer. Then it ' \
            'will map any mated reads against to those refs and also map the nonpaired reads against that ref in a separate ' \
            'call. Then when both are finished it will convert to bam/sort/index/merge/reindex the results. It attempts all of ' \
            'this inside of the /dev/shm filesystem which should be very fast. If /dev/shm cannot be used then /tmp will be used. '\
            'If you want the temporary files that are created to stay then you can use the --keep-temp argument'
    )
    
    parser.add_argument(
        dest='reads',
        help='Location to a directory containing fastq read files'
    )

    parser.add_argument(
        dest='reference',
        help='Location of the reference to map to. Should be a fasta file or directory of fasta files that will be compiled together.'
    )

    valid_platforms = ['MiSeq','Sanger']
    parser.add_argument(
        '--platforms',
        dest='platforms',
        nargs='+',
        choices=valid_platforms,
        default=valid_platforms,
        help='List of platforms to include[Default:{}]'.format(valid_platforms)
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        default='bwa_mem.bam',
        help='Where the output bam should be placed[Default: bwa_mem.bam]'
    )

    parser.add_argument(
        '--keep-temp',
        dest='keep_temp',
        action='store_true',
        default=False,
        help='Flag to indicate that you want the temporary files kept instead of removing them which is the default action'
    )

    return parser.parse_args( args )

class InvalidReference(Exception): pass

def bwa_mem( read1, mate=None, ref=None, output='bwa.sai' ):
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
        ref = compile_refs( ref )
    # First, make sure the reference is indexed
    if not index_ref(ref):
        raise InvalidReference("{} cannot be indexed by bwa")

    # Setup BWA Mem
    mem = None
    if mate:
        mem = BWAMem( ref, read1, mate, bwa_path=which_bwa() )
    else:
        mem = BWAMem( ref, read1, bwa_path=which_bwa() )

    ret = mem.run( output )
    print "Ret: " + str(ret)
    if ret != 0:
        return ret
    else:
        return output
