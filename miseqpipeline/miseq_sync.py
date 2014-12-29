"""
The intention of this script is to easily sync a given MiSeq run path from the sequencer into the :doc:`NGSData <ngsdata>` structure

You need to ensure that the run directory for the run you want to sync is available and that you know the path to it.

Accessing MiSeq Data
====================

The MiSeq stores runs on a separate drive on the MiSeq(Probably the D:\ drive). You will need to ensure that this directory is shared or at least the MiSeqOutput directory is shared. See :ref:`create-share-user`

On the computer that you will be running miseq_sync.py from you will need to ensure that the MiSeq share is mounted somewhere on the system. A good practice is to create a folder somewhere called Instruments and then under there create folders for each of your sequencers.

**Example**

    .. code-block:: bash

        mkdir -p /Instruments/MiSeq

Then you can mount the MiSeq shared drive to that folder.
See :ref:`mount-cifs-linux`

Usage
=====

At this time there is very little output from the miseq_sync.py command until it finishes copying data which can take anywhere from 30 minutes to 2 hours depending on data sizes and network congestion. Be patient and scan through the output to look for failures after it finishes.

    .. code-block:: bash

        miseq_sync.py /path/to/Illumina/MiSeqOutput/runname

Verify Samples Synced
---------------------

coming soon...

How it works
============

The miseq instrument reads the SampleSheet.csv and numbers each row listed under the comma separated values section. It assings each row a number starting with 0 which is the sample number.

Under each MiSeq run's raw data directory inside of Data/Intensities/BaseCalls it creates the gzipped(compressed) fastq file for each paired read.

You can read more about the naming structure `here <http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm>`_

Sync Process
------------

For the examples, we will assume that your NGS Data directory is::

    /home/NGSData

and your miseq run name is::

    140101_M00001_0001_000000000-AAAAA

#. First the script syncs the fastq.gz files from inside of <run name>/Data/Intensities/BaseCalls/ to the destination specified when you run the script using the --ngsdata

    This will create::

        /home/NGSData/RawData/MiSeq/140101_M00001_0001_000000000-AAAAA/Data/Intensities/BaseCalls

    and then copies all .fastq.gz from your run into there

    You can view the `FASTQ Files <http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm>`_ site to see how the MiSeq names it's .fastq.gz files

    An example paired .fastq.gz would be::

        SAMPLENAME1_S1_L001_R1_001.fastq.gz
        SAMPLENAME1_S1_L001_R2_001.fastq.gz

    This name can then be parsed as follows:

        * SAMPLENAME1 is the name of the sample
        * S1 indicates it is sample #1(aka first sample in the samplesheet)
        * L001 means it was lane 1
            * Most likely all sample files will have L001
        * R1 means Read 1(forward), R2 means Read 2(reverse)
        * 001 then indicates it was set 1
            * Most likely all sample files will have 001 unless you changed some setting related to ``--fastq-cluster-count``
    
#. The script then creates a directory with the same name as the MiSeq raw data directory under ReadData/MiSeq

    Creates::

        /home/NGSData/ReadData/140101_M00001_0001_000000000-AAAAA

#. Then it unpacks each fastq.gz file into that directory so there are only .fastq files and appends the run date to them before the .fastq

    So in the previous example the two file names would become::

        SAMPLENAME1_S1_L001_R1_001_2014_01_01.fastq
        SAMPLENAME1_S1_L001_R2_001_2014_01_01.fastq

    and be placed inside of::

        /home/NGSData/ReadData/140101_M00001_0001_000000000-AAAAA

#. Once the files are unpacked it uses the <sample name> field in the each of the file names to determine the final samplename to use

        Aka::

            SAMPLENAME1

#. It then ensures there is a directory with each samplename under ReadsBySample

        Per our example, the following directory would be created::
            
            /home/NGSData/ReadsBySample/SAMPLENAME1

#. It then creates a symlink with the same name as the file in the ReadData/MiSeq/<rundirectory> inside of the appropriate ReadsBySample/<samplename> that links back to the .fastq file

    The symlinks would come out as follows::

        SAMPLENAME1_S1_L001_R1_001_2014_01_01.fastq -> ../../ReadData/MiSeq/140101_M00001_0001_000000000-AAAAA/SAMPLE1_S01_L001_R1_001_2014_01_01.fastq
        SAMPLENAME1_S1_L001_R2_001_2014_01_01.fastq -> ../../ReadData/MiSeq/140101_M00001_0001_000000000-AAAAA/SAMPLE1_S01_L001_R2_001_2014_01_01.fastq
        
#. Once all ReadsBySample directories and files are symlinked it finishes syncing the rest of the miseq run directory

File Naming Scheme
==================


"""

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
    # Ensure all can be converted to integers
    try:
        int(y)
        int(m)
        int(d)
    except ValueError as e:
        raise ValueError( 'run directory does not contain a valid date in the beginning' )
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
    # Fixes #1060. Ensures no trailing slash
    src = normpath( src )
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
    from miseqpipeline import config
    conf_parser, args, config, configfile = config.get_config_argparse(args)
    defaults = config['miseq_sync']

    parser = argparse.ArgumentParser(
        description='Sync MiSeq run into the NGSData structure',
        parents=[conf_parser]
    )

    parser.add_argument(
        'runpath',
        help='Path to the run to be synced(/path/to/YYMMDD_)'
    )

    parser.add_argument(
        '--ngsdata',
        default=defaults['ngsdata']['default'],
        help=defaults['ngsdata']['help']
    )

    return parser.parse_args( args )

