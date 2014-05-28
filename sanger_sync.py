import shutil
from os.path import *
import os
from glob import glob
from Bio import SeqIO
import re

import log
logger = log.setup_logger( __name__, log.get_config() )

# For invalid formatted filenames
class InvalidFormat(Exception): pass

def sync_sanger( runpath, ngsdata ):
    pass

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
            logger.info( '{} already existed so it was skipped'.format(rdpath) )
        fqpath = rdpath.replace('.ab1', '.fastq' )
        if not exists( fqpath ):
            logger.info( 'Converting {} to fastq {}'.format(rdpath,fqpath) )
            SeqIO.convert( rdpath, 'abi', fqpath, 'fastq' )
        else:
            logger.info( '{} already existed so it was skipped'.format(fqpath) )

def samplename_from_read( filepath ):
    p = '(\S+?)_[FR]\d+_\d{4}_\d{2}_\d{2}_\S+?_\S+?_[A-H]\d{2}.(fastq|ab1)'
    m = re.match( p, basename( filepath ) )
    if not m:
        raise InvalidFormat( '{} is not a valid Sanger filename'.format(filepath) )
    return m.groups(0)[0]











'''
 for r in reads:
...   sn = re.match( '(\S+?)_[FR]\d*', basename(r) ).groups(1)[0]
...   run = dirname( r )
...   lnk = relpath( r, join('/home','EIDRUdata','NGSData', 'ReadsBySample', sn ) )
...   sd = join( '/home', 'EIDRUdata', 'NGSData', 'ReadsBySample', sn )
...   if not isdir( sd ):
...     os.makedirs( sd )
...   os.symlink( lnk, join( sd, basename(r) ) )

'''
