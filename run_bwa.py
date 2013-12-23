from bwa.bwa import BWAMem, index_ref, which_bwa, compile_refs

import os
import sys
from os.path import *
import logging

log = logging.getLogger(__name__)

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
            'this inside of the /dev/shm filesystem which should be very fast. If you want the temporary files that are created to stay '\
            'then you can use the --keep-temp argument'
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
