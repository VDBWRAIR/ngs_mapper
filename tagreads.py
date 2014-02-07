#!/usr/bin/env python

import argparse
import sys
import pysam
import re
import shutil
import os.path

import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
log = logging.getLogger('tagreads')

# Exception for when headers exist
class HeaderExists(Exception): pass

# The next 3 tuples have to be the same length and each index in each is related the same index in each tuple
# AKA zip( IDS, PLATFORMS, ID_MAP ) should work as expected
# Read group ID list
IDS = ('Roche454', 'IonTorrent', 'MiSeq', 'Sanger')
# Valid platforms for read groups
PLATFORMS = ('L454', 'IONTORRENT', 'ILLUMINA', 'CAPILLARY')
# Read name map to ID name
ID_MAP = (
    re.compile( '[0-9A-Z]{14}' ),
    re.compile( '[A-Z0-9]{5}:\d{5}:\d{5}' ),
    re.compile( 'M[0-9]{5}:\d:\d{9}-[A-Z0-9]{5}:\d:\d{4}:\d{4,5}:\d{4,5}' ),
    re.compile( '.*' )
)
# Read Group Template
RG_TEMPLATE = {
    'SM': None,
    'ID': None,
    'PL': None,
    'CN': None
}

def main( args ):
    for bam in args.bamfiles:
        tag_bam( bam, args.SM, args.CN )

def tag_bam( bam, SM, CN ):
    log.info( "Gathering existing header for {}".format(bam) )
    hdr = get_rg_headers( bam, SM, CN )
    tag_reads( bam, hdr )

def tag_reads( bam, hdr ):
    ''' Sets header of bam and tags all reads appropriately for each platform '''
    # Open the existing bam to fetch the reads to modify from
    untagged_bam = pysam.Samfile( bam )
    # Start a new bam file with new header
    if 'RG' in untagged_bam.header:
        raise HeaderExists( "Refusing to add header as RG header already exists." )
    tagged_bam = pysam.Samfile( bam + '.tagged', 'wb', header = hdr )
    # Tag the reads
    log.info( "Tagging reads for {}".format(bam) )
    for read in untagged_bam.fetch(until_eof=True):
        rg = get_rg_for_read( read )
        # Have to modify it to get it to save according to docs
        read.tags = read.tags + [('RG',rg)]
        tagged_bam.write( read )
    tagged_bam.close()
    untagged_bam.close()
    log.info( "Finished tagging reads for {}".format(bam) )
    log.info( "Sorting {}".format(bam) )
    pysam.sort( bam + '.tagged', bam.replace('.bam','') )
    log.info( "Indexing {}".format(bam) )
    pysam.index( bam )

def get_rg_for_read( aread ):
    ''' Gets the read group name for the given pysam.AlignedRead '''
    rname = aread.qname
    for i, p in enumerate( ID_MAP ):
        if p.match( rname ):
            return IDS[i]
    raise UnknownReadNameFormat( "{} is from an unknown platform and cannot be tagged".format(rname) )

def get_rg_headers( bam, SM=None, CN=None ):
    old_header = get_bam_header( bam )
    old_header['RG'] = []

    for id, pl in zip( IDS, PLATFORMS ):
        rg = RG_TEMPLATE.copy()
        if SM is not None:
            rg['SM'] = SM
        else:
            rg['SM'] = os.path.basename(bam).replace( '.bam', '' )
        if CN is not None: rg['CN'] = CN
        rg['ID'] = id
        rg['PL'] = pl
        old_header['RG'].append( rg )

    return old_header

def get_bam_header( bam ):
    s = pysam.Samfile( bam )
    hdr = s.header
    s.close()
    return hdr

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='''Adds Sanger, MiSeq, 454 and IonTorrent read groups to a bam file.
        Will tag all reads in the alignment for each of the platforms based on the read name'''
    )

    parser.add_argument(
        dest='bamfiles',
        nargs='+',
        help='Bamfile to add read groups to for platforms'
    )

    parser.add_argument(
        '-SM',
        dest='SM',
        default=None,
        help='Sets the SM tag value inside of each read ' \
            'group to the value specified. Default is to use the portion of the filename ' \
            'that precedes the .bam'
    )

    parser.add_argument(
        '-CN',
        dest='CN',
        default=None,
        help='Sets the CN tag value inside of each read ' \
            'group to the value specified. Default is to not include the CN tag'
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main( parse_args() )
