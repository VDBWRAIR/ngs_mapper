"""
The intention of this script is to easily sync a given Sanger run path from the sequencer into the :doc:`NGSData <../ngsdata>` structure

You need to ensure that the run directory for the run you want to sync is available and that you know the path to it.

Accessing Sanger Data
=====================

You will need to configure the Sanger instrument to save your files in a specific format.

You will also have to ensure that the location that your Sanger instrument saves its data to is :ref:`shared <create-share-user>`

On the computer that you will be running sanger_sync from you will need to ensure that the Sanger share is mounted somewhere on the system. A good practice is to create a folder somewhere called Instruments and then under there create folders for each of your sequencers.

**Example**

    .. code-block:: bash

        mkdir -p /Instruments/Sanger

Then you can mount the Sanger shared drive to that folder.
See :ref:`mount-cifs-linux`

Usage
=====

At this time there is very little output from the sanger_sync command until it finishes copying data which can take anywhere from 30 minutes to 2 hours depending on data sizes and network congestion. Be patient and scan through the output to look for failures after it finishes.

    .. code-block:: bash

        sanger_sync /path/to/Sanger/Run_3130xl...

Verify Samples Synced
---------------------

coming soon...

How it works
============

#. Copy Run_3130xl directory(only including .ab1 files) to RawData/Sanger/
#. Create Run_3130xl directory under ReadData/Sanger/ with same name as original Run_3130xl directory
    #. Symlink all original .ab1 files into this directory
    #. Convert all .ab1 to .fastq 
#. Parse the sanger filename and create ReadsBySample/samplename directory
#. Symlink all .fastq and .ab1 files for that samplename from ReadData into Samplename directory

"""
import shutil
from os.path import *
import os
from glob import glob
from Bio import SeqIO
import re
import sys

import log
logger = log.setup_logger( basename(__file__), log.get_config() )

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
            logger.info( 'Copied {0} to {1}'.format(read,rp) )
            shutil.copy( read, rp )
            synced.append( rp )
        else:
            logger.info( '{0} already existed so it was skipped'.format(rp) )
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
            logger.info( 'Symlinking {0} to {1}'.format(rdpath, lnk) )
            cd = os.getcwd()
            os.symlink( lnk, rdpath )
        else:
            logger.info( 'Skipping existing abi file {0}'.format(rdpath) )
        fqpath = rdpath.replace('.ab1', '.fastq' )
        if not exists( fqpath ):
            logger.info( 'Converting {0} to fastq {1}'.format(rdpath,fqpath) )
            SeqIO.convert( rdpath, 'abi', fqpath, 'fastq' )
        else:
            logger.info( 'Skipping existing fastq file {0}'.format(fqpath) )

def samplename_from_read( filepath ):
    p = '(\S+?)_[FR]\d+_\d{4}_\d{2}_\d{2}_\S+?_\S+?_[A-H]\d{2}.(fastq|ab1)'
    m = re.match( p, basename( filepath ) )
    if not m:
        raise InvalidFormat( '{0} is not a valid Sanger filename'.format(filepath) )
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
            logger.info( 'Symlinking {0} to {1}'.format(rdpath, lnk) )
            os.symlink( lnk, rdpath )
        else:
            logger.info( 'Skipping existing file {0}'.format(rdpath) )

def main():
    args = parse_args()
    sync_sanger( args.runpath, args.ngsdata )

def parse_args( args=sys.argv[1:] ):
    import argparse

    from ngs_mapper import config
    conf_parser, args, config, configfile = config.get_config_argparse(args)
    defaults = config['sanger_sync']

    parser = argparse.ArgumentParser(
        description='Syncs Sanger Run_ directories into the NGSData structure',
        parents=[conf_parser]
    )

    
    parser.add_argument(
        '--ngsdata',
        dest='ngsdata',
        default=defaults['ngsdata']['default'],
        help=defaults['ngsdata']['help']
    )

    parser.add_argument(
        'runpath',
        help='Path to Sanger Run_3130xl directory'
    )

    return parser.parse_args( args )
