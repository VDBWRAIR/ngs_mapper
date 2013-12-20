from bwa import seqio
from bwa.bwa import BWAMem, index_ref, which_bwa
#import bwa.bwa

import os
from os.path import *
import logging

log = logging.getLogger(__name__)

class InvalidReadFile(Exception): pass

def compile_reads( readfilelist, outputdir ):
    '''
        Compiles all read files inside of readfilelist into respective files.
        Creates F.fq, R.fq and/or NP.fq depending on the reads found in readfilelist
        Only compiles fastq files. If others are given an exception will be raised

        @param readfilelist - List of read file paths. If any of the items are a 2 item tuple, that item is treated as a
            mate pair read set and the first item will be the Forward read file and the second the Reverse read file.
        @param outputdir - Where to output the three files

        @returns a dictionary {'F': join(outputdir,'F.fq'), 'R': join(outputdir,'R.fq'), 'NP': join(outputdir,'NP.fq')}
            If there were no mated files given then F & R will be None. Same goes for NP
    '''
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    files_written = {'F':[],'R':[],'NP':[]}

    # Build the list of files to be written for each type
    for read in readfilelist:
        # Non Paired read
        if isinstance(read,str):
            if read.endswith('.fastq'):
                files_written['NP'].append(read)
            else:
                raise InvalidReadFile("{} is not a fastq file. Only fastq files are supported at this time.".format(read))
        elif isinstance(read,tuple) and len(read) == 2:
            if read[0].endswith('.fastq'):
                files_written['F'].append(read[0])
            else:
                raise InvalidReadFile("{} is not a fastq file. Only fastq files are supported at this time.".format(read[0]))
            if read[1].endswith('.fastq'):
                files_written['R'].append(read[1])
            else:
                raise InvalidReadFile("{} is not a fastq file. Only fastq files are supported at this time.".format(read[1]))
        else:
            raise ValueError("Somehow neither got 1 or 2 items for a read in readfilelist")
    
    # Now concat the files to their respective output file
    for f,files in files_written.items():
        if files:
            outfile = join(outputdir,f+'.fq')
            seqio.concat_files(files, outfile)
            files_written[f] = outfile
        else:
            files_written[f] = None

    return files_written

class InvalidReference(Exception): pass

def bwa_mem( read1, mate=None, ref=None, output='bwa.sai' ):
    '''
        Runs the bwa mem algorithm on read1 against ref. If mate is given then run that file with the read1 file
        so paired alignment is done.

        TODO:
            bwa_path should be an option to specify where the executable is

        @param read1 - File path to read
        @param mate - Mate file path
        @param ref - Reference file path
        @param output - The output destination

        @returns the output path if sucessful or -1 if something went wrong
    '''
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
