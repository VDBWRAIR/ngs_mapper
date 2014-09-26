import argparse
import os
import sys
from os.path import *
import subprocess
import shlex
import shutil
from glob import glob
from datetime import datetime
import gzip
import csv

import log
logger = log.setup_logger( basename(__file__), log.get_config() )

def main( args ):
    sync( args.runpath, args.ngsdata )

#Tested
def file_already_copied( src, dst ):
    '''
        Just check if stat.st_size are same
    '''
    return os.stat( src ).st_size == os.stat( dst ).st_size

#Tested
def get_basecalls_dir( runpath ):
    '''
        Get path to BaseCalls
    '''
    # BaseCalls relative path(contains fastq.gz)
    return join( runpath, 'Data', 'Intensities', 'BaseCalls' )

#Tested
def get_rundate( rundir ):
    '''
        Get the date for a rundir by extracting the first 6 digits
    '''
    rundir = basename( rundir )
    y = rundir[0:2]
    m = rundir[2:4]
    d = rundir[4:6]
    return '20{}_{}_{}'.format(y,m,d)

#Tested
def samplename_from_fq( fastqp ):
    '''
        Return sample name from miseq fastq path
    '''
    fastqp = basename( fastqp )
    return fastqp.split( '_' )[0]

def parse_samplesheet( sheetpath ):
    '''
        Parses sample sheet returning list of everything under [Data] as a csv.DictReader
    '''
    with open( sheetpath ) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith( '[Data]' ):
                # Eat the header and get to the good stuff
                break
        for sample in csv.DictReader( fh ):
            yield sample

def sync_fastq( srcrun, ngsdata ):
    src_run_path = srcrun
    runname = basename( srcrun )
    src_fastq_path = get_basecalls_dir( src_run_path )
    dst_run_path = join( ngsdata, 'RawData', 'MiSeq', runname )
    dst_fastq_path = get_basecalls_dir( dst_run_path )
    if not exists( dst_fastq_path ):
        os.makedirs( dst_fastq_path )
    gz = glob( join( src_fastq_path, '*.fastq.gz' ) )
    for fq in gz:
        dst = join( dst_fastq_path, basename( fq ) )
        if exists( dst ) and file_already_copied( fq, dst ):
            logger.info( '{} already exists at {}'.format(fq, dst_fastq_path) )
        else:
            logger.info( 'Copying {} into {}'.format(fq, dst_fastq_path) )
            shutil.copy( fq, dst_fastq_path )

def sync( src, ngsdata ):
    '''
        Sync all files/dirs from src into ngsdata
    '''
    sync_fastq( src, ngsdata )
    srundir = join( ngsdata, 'RawData', 'MiSeq', basename(src) )
    create_readdata( srundir, ngsdata )
    link_reads( srundir, ngsdata )
    rsync_run( src, ngsdata )

def rsync_run( rundir, ngsdata ):
    # No trailing /
    src = normpath( rundir )
    # Dst will end in MiSeq since we removed the trailing / in src
    dst = join( ngsdata, 'RawData', 'MiSeq' )
    logger.info( 'The fastq read data is synced. Syncing the rest of the data' )
    cmd = 'rsync -av --progress --size-only {} {}'.format( src, dst )
    cmd = shlex.split( cmd )
    logger.debug( subprocess.check_output( cmd ) )


def create_readdata( rundir, ngsdata ):
    '''
        Uncompress fastq files from rundir into basename(rundir)'s ReadData
        and add date to the end
    '''
    fqpath = get_basecalls_dir( rundir )
    dstroot = join( ngsdata, 'ReadData', 'MiSeq', basename( rundir ) )
    rundate = get_rundate( rundir )
    if not exists( dstroot ):
        os.makedirs( dstroot )
    for gz in glob( join( fqpath, '*.fastq.gz' ) ):
        dstfq = join( dstroot, basename( gz ).replace( '.fastq.gz', '_{}.fastq'.format( rundate ) ) )
        if not exists( dstfq ):
            logger.info( 'Unpacking {} to {}'.format(gz, dstfq) )
            with gzip.open( gz, 'rb' ) as fr, open( dstfq, 'w' ) as fw:
                fw.write( fr.read() )
        else:
            logger.debug( '{} looks to be unpacked already as {}'.format(gz, dstfq) )

def link_reads( rundir, ngsdata ):
    '''
        Symlink read files into ReadsBySample
    '''
    readdata = join( ngsdata, 'ReadData', 'MiSeq', basename( rundir ) )
    readsbysample = join( ngsdata, 'ReadsBySample' )
    for fq in glob( join( readdata, '*.fastq' ) ):
        samplename = samplename_from_fq( fq )
        # Ensure sampledir
        sampledir = join( readsbysample, samplename )
        if not exists( sampledir ):
            os.makedirs( sampledir )
        rpath = relpath( fq, sampledir )
        lnkdst = join( sampledir, basename( fq ) )
        if not exists( lnkdst ):
            logger.info( 'Symlinking {} to {}'.format( rpath, lnkdst ) )
            os.symlink( rpath, lnkdst )
        else:
            logger.debug( '{} already exists.'.format( lnkdst ) )

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Sync MiSeq run into the NGSData structure'
    )

    parser.add_argument(
        'runpath',
        help='Path to the run to be synced(/path/to/YYMMDD_)'
    )

    ngsdata_default = '/home/EIDRUdata/NGSData'
    parser.add_argument(
        '--ngsdata',
        default=ngsdata_default,
        help='Path to the root of the NGSData data structure[Default: {}]'.format(ngsdata_default)
    )

    return parser.parse_args( args )

