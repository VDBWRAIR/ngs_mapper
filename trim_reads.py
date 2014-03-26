#!/usr/bin/env python

import subprocess
import os
import argparse
import sys
from os.path import basename, join, isdir
from glob import glob

import log
lconfig = log.get_config()
logger = log.setup_logger( 'trim_reads', lconfig )

def main( args ):
    trim_reads_in_dir(
        args.readsdir,
        args.q,
        args.outputdir
    )

def trim_reads_in_dir( readdir, qual_th, out_path ):
    '''
        Trims all read files in a given directory and places the resulting files into out_path directory

        @param readdir - Directory with read files in it(sff and fastq only)
        @qual_th - What to pass to cutadapt -q
        @out_path - Output directory path
    '''
    # Only sff and fastq files
    reads = [f for f in os.listdir(readdir) if f.endswith('sff') or f.endswith('fastq')]
    # Make out_path
    if not isdir( out_path ):
        os.mkdir( out_path )
    # Trim all the reads
    for read in reads:
        reado = read.replace('.sff','.fastq')
        trim_read( join(readdir,read), qual_th, join(out_path,reado) )

def trim_read( readpath, qual_th, out_path=None ):
    '''
        Trims the given readpath file and places it in out_path
        If out_path not given then just put it in current directory with the same basename

        @param readpath - Path to the read to trim .fastq and .sff support only
        @param qual_th - Quality threshold to trim reads on
        @param out_path - Where to put the trimmed file

        @returns path to the trimmed fastq file
    '''
    from Bio import SeqIO
    tfile = None
    if out_path is None:
        out_path = basename( readpath ).replace('.sff','.fastq')
    logger.debug( "Using {} as the output path".format(out_path) )

    # Convert sff to fastq
    if readpath.endswith('.sff'):
        logger.debug( "Converting {} to fastq".format(readpath) )
        # Just put in temp location then remove later
        tfile = '/tmp/sff.fastq'
        try:
            SeqIO.convert( readpath, 'sff', tfile, 'fastq' )
        except AssertionError as e:
            # Ignore the stupid sff bug in BioPython
            pass
        readpath = tfile

    # Run cutadapt on the file
    stats = run_cutadapt( readpath, out_path, q=qual_th )

    # Clean up temp file
    if tfile:
        os.unlink(tfile)

    return out_path

def run_cutadapt( *args, **kwargs ):
    '''
        Runs cutadapt with the given arguments and kwargs
        
        @param - fastq file to trim
        @param - output file location
        @param q - Quality threshold

        @returns the stderr output from cutadapt
    '''
    cmd = ['cutadapt', '-q', str(kwargs.get('q')), args[0]]
    fout = args[1]
    if isinstance( fout, str ):
        fout = open(fout,'wb')
    # Write stdout to output argument(should be fastq)
    # Allow us to read stderr which should be stats from cutadapt
    logger.debug( "Running {}".format(cmd) )
    p = subprocess.Popen( cmd, stdout=fout, stderr=subprocess.PIPE )
    # Only stderr should be available
    _,se = p.communicate()
    if p.returncode != 0:
        e = subprocess.CalledProcessError(p.returncode,' '.join(cmd), se)
        raise e
    return se

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Trims reads'
    )

    parser.add_argument(
        dest='readsdir',
        help='Read or directory of read files'
    )

    qual_default=20
    parser.add_argument(
        '-q',
        dest='q',
        default=qual_default,
        help='Quality threshold to trim[Default:{}]'.format(qual_default)
    )

    outputdir_default='trimmed_reads'
    parser.add_argument(
        '-o',
        dest='outputdir',
        default=outputdir_default,
        help='Where to output the resulting files[Default:{}]'.format(outputdir_default)
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main(parse_args())
