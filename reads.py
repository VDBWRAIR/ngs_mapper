import os
from os.path import *

import logging
log = logging.getLogger(__name__)

# Valid Read extensions
VALID_READ_EXT = ('sff','fastq')

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
    from bwa import seqio
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    files_written = {'F':[],'R':[],'NP':[]}

    # Build the list of files to be written for each type
    for read in readfilelist:
        # Non Paired read
        if isinstance(read,str):
            if is_valid_read( read ):
                files_written['NP'].append(read)
            else:
                raise InvalidReadFile("{} is not a fastq file. Only fastq files are supported at this time.".format(read))
        elif isinstance(read,tuple) and len(read) == 2:
            if is_valid_read( read[0] ):
                files_written['F'].append(read[0])
            else:
                raise InvalidReadFile("{} is not a fastq file. Only fastq files are supported at this time.".format(read[0]))
            if is_valid_read( read[1] ):
                files_written['R'].append(read[1])
            else:
                raise InvalidReadFile("{} is not a fastq file. Only fastq files are supported at this time.".format(read[1]))
        else:
            raise ValueError("Somehow neither got 1 or 2 items for a read in readfilelist")

    # Now concat the files to their respective output file
    for f,files in files_written.items():
        if files:
            outfile = join(outputdir,f+'.fq')
            log.info( "Compiling reads from {} into {}".format( files, outfile ) )
            # List of only sff
            sffs = [r for r in files if r.endswith('.sff')]
            # List of all others
            other = [r for r in files if not r.endswith('.sff')]

            # All files to concat together
            cfiles = []

            # Cat all sff
            if sffs:
                sff_out = join(outputdir,'sffs.fq')
                log.debug( "Converting {} to temp fastq file {}".format(sffs,sff_out) )
                try:
                    r = seqio.sffs_to_fastq( sffs, sff_out )
                except AssertionError as e:
                    # Skip this error that started occuring in new version of BioPython
                    # Pretty big assumption that the outfiles are ok
                    # Check Biopython issue #294
                    pass
                cfiles.append( sff_out )
            if other:
                # Just concat these as normal since they are just text
                cfiles += other

            # concat everything together now
            log.debug( "Concating all files {} to {}".format(cfiles,outfile) )
            seqio.concat_files(cfiles, outfile)
            # Remove the sff concatted file
            if sffs:
                os.unlink( sff_out )

            # Record the file that was written
            files_written[f] = outfile
        else:
            files_written[f] = None

    return files_written

def is_valid_read( readpath ):
    '''
        Just checks to make sure file extension is in VALID_READ_EXT

        @param readpath - Path to read file

        @returns true if ext of readpath is in VALID_READ_EXT
    '''
    # File extension without the period in VALID_READ_EXT ?
    return splitext( readpath )[1][1:] in VALID_READ_EXT
