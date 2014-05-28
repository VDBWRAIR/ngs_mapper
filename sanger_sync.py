#!/usr/bin/env python

import shutil
from os.path import *
import os
from glob import glob
from Bio import SeqIO
import re
import sys

import log
logger = log.setup_logger( __name__, log.get_config() )

# For invalid formatted filenames
class InvalidFormat(Exception): pass

def sync_sanger( runpath, ngsdata ):
    rund = basename( runpath )
    rawd = join( ngsdata, 'RawData', 'Sanger', rund )
    readd = join( ngsdata, 'ReadData', 'Sanger', rund )

    sync_run( runpath, ngsdata )
    sync_readdata( rawd, ngsdata )
    link_reads( readd, ngsdata )

def sync_run( runpath, ngsdata ):
    '''
    Ensures that all ab1 files are copied into the RawData/Sanger
    directory under the Run directory's folder and that filenames are correct

    @param runpath - Path to the run to be synced
    @param ngsdata - Root NGSData path

    @returns list of full paths to all synced files(only files that were transfered)
    '''
    rawd = join( ngsdata, 'RawData', 'Sanger' )
    run = basename(runpath)
    rund = join( rawd, run )

    if not isdir( rund ):
        os.makedirs( rund )

    reads_to_copy = glob( join( runpath, '*.ab1' ) )
    # Should raise error if any of the filenames have incorrect formats
    [samplename_from_read(r) for r in reads_to_copy]
    synced = []
    for read in reads_to_copy:
        rp = join( rund, basename(read) )
        if not exists( rp ):
            logger.info( 'Copied {} to {}'.format(read,rp) )
            shutil.copy( read, rp )
            synced.append( rp )
        else:
            logger.info( '{} already existed so it was skipped'.format(rp) )
    return synced

def sync_readdata( rawdir, ngsdata ):
    '''
    Ensures that ab1 files are symlinked from the RawData/Sanger/Run directory
    and that they are then converted to fastq
    
    @param rawdir - RawData/Sanger/Run path
    @param ngsdata - Path to root NGSData directory
    '''
    raw_reads = glob( join( rawdir, '*.ab1' ) )
    readd = join( ngsdata, 'ReadData', 'Sanger', basename(rawdir) )
    if not isdir( readd ):
        os.makedirs( readd )
    for read in raw_reads:
        lnk = relpath( read, readd )
        rdpath = join( readd, basename(read) )
        if not exists( rdpath ):
            logger.info( 'Symlinking {} to {}'.format(rdpath, lnk) )
            cd = os.getcwd()
            os.symlink( lnk, rdpath )
        else:
            logger.info( 'Skipping existing abi file {}'.format(rdpath) )
        fqpath = rdpath.replace('.ab1', '.fastq' )
        if not exists( fqpath ):
            logger.info( 'Converting {} to fastq {}'.format(rdpath,fqpath) )
            SeqIO.convert( rdpath, 'abi', fqpath, 'fastq' )
        else:
            logger.info( 'Skipping existing fastq file {}'.format(fqpath) )

def samplename_from_read( filepath ):
    p = '(\S+?)_[FR]\d+_\d{4}_\d{2}_\d{2}_\S+?_\S+?_[A-H]\d{2}.(fastq|ab1)'
    m = re.match( p, basename( filepath ) )
    if not m:
        raise InvalidFormat( '{} is not a valid Sanger filename'.format(filepath) )
    return m.groups(0)[0]

def link_reads( readdata, ngsdata ):
    '''
    Ensures that all files from readdata are symlinked into ReadsBySample/SampleName/
    Does not overwrite existing links

    @param readdata - ReadData/Sanger/Rund path
    @param ngsdata - Root path to NGSData directory
    '''
    read_files = glob( join( readdata, '*' ) )
    rbs_root = join( ngsdata, 'ReadsBySample' )
    for read in read_files:
        sn = samplename_from_read( read )
        rbs = join( rbs_root, sn )
        lnk = relpath( read, rbs )
        rdpath = join( rbs, basename( read ) )
        if not isdir( rbs ):
            os.makedirs( rbs )
        if not islink( rdpath ):
            logger.info( 'Symlinking {} to {}'.format(rdpath, lnk) )
            os.symlink( lnk, rdpath )
        else:
            logger.info( 'Skipping existing file {}'.format(rdpath) )

def main( args ):
    sync_sanger( args.runpath, args.ngsdata )

def parse_args( args=sys.argv[1:] ):
    import argparse
    parser = argparse.ArgumentParser(
        description='Syncs Sanger Run_ directories into the NGSData structure'
    )

    ngsdata = '/home/EIDRUdata/NGSData'
    parser.add_argument(
        '--ngsdata',
        dest='ngsdata',
        default=ngsdata,
        help='NGSData root path'
    )

    parser.add_argument(
        'runpath',
        help='Path to Sanger Run_3130xl directory'
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main( parse_args() )
