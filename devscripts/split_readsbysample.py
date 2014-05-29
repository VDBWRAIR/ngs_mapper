import csv
from os.path import *
import os
import sys
import re
from glob import glob

def main( args ):
    # Ensure new RBS is created
    if not isdir( args.outrbs ):
        os.makedirs( args.outrbs )
    index = index_samples( args.samplesheet )
    split_samples( index, args.inrbs, args.outrbs )
    
def split_samples( index, inrbs, outrbs ):
    for samplename, samples in index.iteritems():
        split_sample( samplename, samples, inrbs, outrbs )

def split_sample( samplename, samples, inrbs, outrbs ):
    ''' Split an existing ReadsBySample/sample into all Sample_ID's for it '''
    # Get the sample_ID's for this sample
    for sampleid, samplenum in samples:
        irbs = join( inrbs, samplename )
        orbs = join( outrbs, sampleid )
        # Make sure dest directory exists
        if not isdir( orbs ):
            os.makedirs( orbs )
        else:
            sys.stderr.write( "Skipping {} because {} already exists\n".format(sampleid,orbs) )
            continue
        glb = join( irbs, '{}_S{}_*'.format( samplename, samplenum ) )
        reads = glob( glb )
        if not reads:
            sys.stderr.write( "There are no existing reads in {}\n".format(irbs) )
        for read in reads:
            outpath = join( orbs, basename(read).replace(samplename,sampleid) )
            os.symlink( read, outpath )

def parse_args( args=sys.argv[1:] ):
    import argparse

    parser = argparse.ArgumentParser(
        description='''Creates ReadsBySample directory from a MiSeq SampleSheet.csv
            by splitting an existing ReadsBySample directory set for a MiSeq run into
            directories by the Sample_ID field'''
    )

    parser.add_argument(
        'samplesheet',
        help='Location of SampleSheet.csv'
    )

    inrbs = '/home/EIDRUdata/NGSData/ReadsBySample'
    parser.add_argument(
        '--readsbysample',
        dest='inrbs',
        default=inrbs,
        help='Path to existing ReadsBySample[Default: {}]'.format(inrbs)
    )

    outrbs = 'ReadsBySample'
    parser.add_argument(
        '--outrbs',
        dest='outrbs',
        default=outrbs,
        help='Path to new ReadsBySample directory to create[Default: {}]'.format(outrbs)
    )

    return parser.parse_args( args )

def index_samples( samplesheet ):
    '''
        >>> index = index_samples( 'SampleSheet.csv' )
        >>> s = index['F071']
        >>> assert s[0] == ('F071-F',1)
        >>> assert s[1] == ('F071-M',41)
    '''
    # Index samples by sample name
    sample_index = {}
    with open(samplesheet) as fh:
        for i, row in enumerate( csv.DictReader(fh), 1 ):
            sampleid = row['Sample_ID']
            samplename = row['Sample_Name']
            samplenum = i
            if samplename not in sample_index:
                sample_index[samplename] = []
            sample_index[samplename].append( (sampleid, samplenum) )
    return sample_index

if __name__ == '__main__':
    main( parse_args() )
