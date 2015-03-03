"""
The intention of this script is to easily sync a given Roche 454 run path from the sequencer into the :doc:`NGSData <../ngsdata>` structure

You need to ensure that the run directory for the run you want to sync is available and that you know the path to it.

You need to have the Data Analysis

Accessing Roche 454 Data
========================

There are quite a few different setups for this instrument, so no specific instructions can be created.

Essentially you will want to mount the /data directory on the Roche analysis workstation to your /Instruments/Roche454 directory.

The goal is to be able to browse to the R\_ directory from the workstation that you will be running roche_sync from

Usage
=====

At this time there is very little output from the roche_sync command until it finishes copying data which can take anywhere from 30 minutes to 2 hours depending on data sizes and network congestion. Be patient and scan through the output to look for failures after it finishes.

    .. code-block:: bash

        roche_sync $(pwd) --ngsdata /path/to/NGSData --midparse ~/ngs_mapper/MidParse.conf

Verify Samples Synced
---------------------

coming soon...

How it works
============

Since the Roche instrument does not automatically create a sample to barcode/region mapping file such as the MiSeq's SampleSheet.csv, you will have to manually create one and place it into the top level R\_ directory.

.. _roche-sample-sheet:

Roche Sample Sheet
------------------

The roche sample sheet is simply a mapping of sample name to the region, barcode and primer fasta file used.
The file must be saved as SampleSheet.csv and placed in the R\_ directory of the run

The format is as follows::

    SampleName,Region,Barcode,PrimerFileLocation
    Sample1,1,RL1,/path/to/primer.fasta
    Sample2,2,TI1,/path/to/primer2.fasta

Sync Process
------------

For the examples, we will assume that your NGS Data directory is::

    /home/NGSData

and your roche run name is::

    R_2014_01_01_01_01_01_FLX00000001_Administrator_SOMEPROJECT

and contains the following signalProcessing directory

    D_2014_01_01_02_02_02_FLX00000001_signalProcessing

#. First the script copies the entire R\_ directory to the NGSData structure's RawData

    This will create::

        /home/NGSData/RawData/Roche454/R_2014_01_01_01_01_01_FLX00000001_Administrator_SOMEPROJECT

#. The script then symlinks the D_*signalProcessing directory into your NGSData/ReadData
#. The script then demultiplexes the signalProcessing sff files
    #. Ensure demultiplexed directory exists inside of signalProcessing
    #. For each sample/barcode listed in the created SampleSheet.csv it will demultiplex that sample out of the main .sff files
        into the demultiplexed directory
    
        For example::

            Sample1__1__RL1__2001_01_01__Unk.sff
            Sample2__2__TI1__2001_01_01__Unk.sff

#. The script then symlinks all demultiplexed reads that are in the RawData/R\_/D\_ directory into the ReadsBySample under the Sample Name directory

    .. code-block:: bash

        Sample1/
            Sample1__1__RL1__2014_01_01__Unk.sff -> ../../RawData/Roche454/R_2014_01_01_01_01_01_FLX00000001_Administrator_SOMEPROJECT/D_2014_01_01_02_02_02_FLX00000001_signalProcessing/demultiplexed/Sample1__1__RL1__2014_01_01__Unk.sff

"""
import os
from os.path import *
import sys
from glob import glob
import subprocess
import re
import csv
import shutil

def main():
    args = parse_args()
    dst = join( args.ngsdata, 'RawData', 'Roche454', basename(args.runpath) )
    sync( args.runpath, args.ngsdata )
    symlink_sigproc( dst, args.ngsdata )
    demultiplex_run( dst, args.midparse )
    link_reads( dst, args.ngsdata )

def sync( src, ngsdata ):
    '''
    Sync all files/dirs in src to dst
    '''
    dst = join( ngsdata, 'RawData', 'Roche454', basename( src ) )
    if isdir( dst ):
        print "{} already exists so skipping data sync. If this directory was only " \
            "partially synced you may want to remove everything and start over".format(dst)
        return
    shutil.copytree( src, dst )

def demultiplex_run( rdir, midparse ):
    '''
    Demultiplex run at rdir into outpath
    '''
    sigprocdir = get_sigprocdir( rdir )
    demuldir = join( sigprocdir, 'demultiplexed' )
    sffdir = join( sigprocdir, 'sff' )
    if not isdir( demuldir ):
        os.mkdir( demuldir )
    region_sff = sff_region_map( sffdir )
    for sample in parse_samplesheet( join(rdir,'SampleSheet.csv') ):
        sname = sample['SampleName']
        region = sample['Region']
        barcode = sample['Barcode']
        primer = sample['PrimerFileLocation']
        sff = join( sffdir, region_sff[region] )
        rdate = get_rundate( rdir )
        sample['date'] = rdate
        outread = format_read_name( **sample )
        outread = join( demuldir, outread )
        # Don't recreate
        if not exists( outread ):
            demultiplex_read( sff, outread, barcode, midparse )

def link_reads( rdir, ngsdata ):
    sigprocdir = get_sigprocdir( rdir )
    demuldir = join( sigprocdir, 'demultiplexed' )
    sffs = glob( join( demuldir, '*.sff' ) )
    readsbysample = join( ngsdata, 'ReadsBySample' )
    for sample in parse_samplesheet( join(rdir,'SampleSheet.csv') ):
        samplename = sample['SampleName']
        sampledir = join( readsbysample, samplename )
        sffs = glob( join(demuldir, samplename + '*.sff') )
        # Ensure sample dir
        if not isdir( sampledir ):
            os.makedirs( sampledir )
        for f in sffs:
            s = relpath( f, sampledir )
            dst = join( sampledir, basename(f) )
            if not exists( dst ):
                os.symlink( s, dst )

def symlink_sigproc( rdir, ngsdata ):
    sigprocdir = get_sigprocdir( rdir )
    readdata = join( ngsdata, 'ReadData', 'Roche454' )
    # Ensure ReadData
    try:
        os.makedirs( readdata )
    except OSError as e:
        pass
    s = relpath( sigprocdir, readdata )
    dst = join( readdata,  basename(sigprocdir) )
    if not exists( dst ):
        os.symlink( s, dst )

def get_rundate( rochedir ):
    '''
    Get date from signalProcessing, imageProcessing or R\_ directory
    '''
    return re.search( 'R_(\d{4}_\d{2}_\d{2})', rochedir ).group(1)

def format_read_name( **kwargs ):
    '''samplename, region, barcode, date, pathogen'''
    nfmt = '{SampleName}__{Region}__{Barcode}__{date}__Unk.sff'
    return nfmt.format(**kwargs)

def demultiplex_read( inputsff, outputsff, barcodename, midparse ):
    '''
    Runs something like

    .. code-block:: bash

        sfffile -o {outputsff} -mcf {midparse} {barcodename}@{inputsff}
    '''
    cmd = ['sfffile','-o',outputsff,'-mcf',midparse,'{}@{}'.format(barcodename,inputsff)]
    print ' '.join( cmd )
    subprocess.call( cmd )

def sff_region_map( sffdir ):
    '''
    Get a mapping of region: basename(sfffilename)
    '''
    sffs = glob( join(sffdir, '*.sff') )
    map = {}
    p = re.compile( '\w+?0(\d).sff' )
    for sff in sffs:
        try:
            r = p.search( sff ).group(1)
        except AttributeError as e:
            raise ValueError( "{} is not a standard multiplexed sff file name".format(sff) )
        map[r] = basename(sff)
    return map

def get_sigprocdir( rdir ):
    '''
    Quick way to get signalProcessing directory from inside of R\_ dir
    '''
    try:
        return glob( join(rdir, '*signalProcessing') )[0]
    except IndexError as e:
        raise ValueError( "{} does not contain a signalProcessing directory".format(rdir) )

def parse_samplesheet( sheetpath ):
    '''
    Parse samplesheet
    @param sheetpath - Path to samplesheet
    @returns csv.DictReader iterator
    '''
    if not exists( sheetpath ):
        raise IOError( "{} does not exist".format(sheetpath) )
    with open( sheetpath ) as fh:
        for row in csv.DictReader( fh ):
            yield row

def parse_args( args=sys.argv[1:] ):
    import argparse
    from ngs_mapper import config
    conf_parser, args, config, configfile = config.get_config_argparse(args)
    defaults = config['roche_sync']

    parser = argparse.ArgumentParser(
        description='Sync roche run directories into the NGSData data structure',
        parents=[conf_parser]
    )

    parser.add_argument(
        'runpath',
        help='Path to the Roche run directory(/path/to/R_something)'
    )

    parser.add_argument(
        '--ngsdata',
        default=defaults['ngsdata']['default'],
        help=defaults['ngsdata']['help']
    )

    import pkg_resources
    default = defaults['midparse']['default']
    if default is None:
        default = pkg_resources.resource_filename(__name__, 'MidParse.conf')
    parser.add_argument(
        '--midparse',
        default=default,
        help=defaults['midparse']['help']
    )

    return parser.parse_args( args )
