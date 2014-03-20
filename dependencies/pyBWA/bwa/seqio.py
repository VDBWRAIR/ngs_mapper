from Bio import SeqIO

import os
import sys
import os.path
import glob
import shutil

class EmptyFileError( Exception ):
    pass

def sffs_to_fastq( sffs, output='sff.fastq' ):
    '''
        Given a list of sffs, concat them into a single fastq
        Nothing created if empty list given

        @raises ValueError if invalid sff file is encountered or output is invalid path
        @param sffs - List of sff file paths to convert and concat into fastq
        @param output - Output fastq file path[Default: sff.fastq]
        @return Path to fastq file created or None if empty list given
    '''
    # Has to be a list
    if not isinstance( sffs, list ):
        raise ValueError( "{} is not a list".format(sffs) )

    # Empty list gets ignored
    if not sffs:
        return

    try:
        # Concat all sequences to output file
        with open( output, 'w' ) as fh:
            for sff in sffs:
                    SeqIO.write( SeqIO.parse( sff, 'sff' ), fh, 'fastq' )
    except (OSError, IOError) as e:
        raise ValueError( "{} is not a valid output file".format(output) )
    except ValueError as e:
        raise ValueError( "{} is not a valid sff file".format(sff) )
    
    return output

def get_reads( dir_path ):
    '''
        Return a list of sff and fastq files in a given dir_path

        @raises ValueError if invalid path given
        @param dir_path - Path to directory of fastq and sff files
        @return list of fastq and sff files found in dir_path. Each has dir_path prefixed to them. Empty list if none found
    '''
    if not os.path.isdir( dir_path ):
        raise ValueError( "{} is not a valid directory".format(dir_path) )
    return glob.glob( os.path.join( dir_path, '*.sff' ) ) + glob.glob( os.path.join( dir_path, '*.fastq' ) ) 

def concat_files( filelist, outputfile ):
    '''
        Duplicate cat *filelist > outputfile
        Don't forget that the files in filelist could end with 2 newlines and thus put empty lines
        into your concatted file. Could be painful with fasta, fastq files
        
        @raises OSError if any fo filelist or outputfile cannot be read/written
        @raises EmptyFileError if outputfile ends up empty
        @raises ValueError if any of filelist are not valid files or if outputfile is not valid

        @param filelist - List of valid files to cat into outputfile.
        @param outputfile - File path to put concatted output into
    '''
    if not isinstance( filelist, list ) or len( filelist ) == 0:
        raise ValueError( "{} is not a valid list of files to concat".format(filelist) )

    if not isinstance( outputfile, str ):
        raise ValueError( "{} is not a valid output path".format(outputfile) )

    if outputfile in filelist:
        raise ValueError( "{} contains the outputfile".format(filelist) )

    # Concat all the found files
    # Could raise IOError or OSError as we are using open on files
    # that are not checked to see if they have perms to read/write
    with open( outputfile, 'wb' ) as fh:
        for f in filelist:
            try:
                with open( f, 'rb' ) as fr:
                    shutil.copyfileobj( fr, fh )
            except (IOError,OSError) as e:
                if e.errno == 2:
                    raise ValueError( "{} does not exist".format(f) )
                os.unlink( outputfile )
                raise e

    if os.stat( outputfile ).st_size == 0:
        os.unlink( outputfile )
        raise EmptyFileError( "Empty files given to concat" )

def seqfile_type( filename ):
    ftype = 'fasta'
    with open( filename ) as fh:
        firstline = fh.readline()
        if firstline.startswith( '>' ):
            ftype = 'fasta'
        elif firstline.startswith( '@' ):
            ftype = 'fastq'
        elif firstline.startswith( '.sff' ):
            ftype = 'sff'
        else:
            raise ValueError( "{} not a valid sequence file".format(
                    filename
                )
            )
    return ftype

def reads_in_file( filename ):
    ftype = seqfile_type( filename )
    return sum( [1 for seq in SeqIO.parse( filename, ftype )] )

