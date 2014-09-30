import os
from os.path import *
import sys
from glob import glob
import subprocess
import re
import csv
import shutil

def main( args ):
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
        Get date from signalProcessing, imageProcessing or R_ directory
    '''
    return re.search( 'R_(\d{4}_\d{2}_\d{2})', rochedir ).group(1)

def format_read_name( **kwargs ):
    '''samplename, region, barcode, date, pathogen'''
    nfmt = '{SampleName}__{Region}__{Barcode}__{date}__Unk.sff'
    return nfmt.format(**kwargs)

def demultiplex_read( inputsff, outputsff, barcodename, midparse ):
    '''
        Runs something like
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
        Quick way to get signalProcessing directory from inside of R_ dir
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
    parser = argparse.ArgumentParser(
        description='Sync roche run directories into the NGSData data structure'
    )

    parser.add_argument(
        'runpath',
        help='Path to the Roche run directory(/path/to/R_something)'
    )

    ngsdata_default = '/home/EIDRUdata/NGSData'
    parser.add_argument(
        '--ngsdata',
        default=ngsdata_default,
        help='Path to the root of the NGSData data structure[Default: {}]'.format(ngsdata_default)
    )

    midparse_default = 'MidParse.conf'
    parser.add_argument(
        '--midparse',
        default=midparse_default,
        help='Path to MidParse.conf file[Default: {}]'.format(midparse_default)
    )

    return parser.parse_args( args )
