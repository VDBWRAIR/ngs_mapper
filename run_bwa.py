from bwa.bwa import BWAMem, index_ref, which_bwa, compile_refs

import os
from os.path import *
import logging

log = logging.getLogger(__name__)

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
