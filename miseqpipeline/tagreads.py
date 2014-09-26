import argparse
import sys
import samtools
import re
import shutil
import os.path
from miseqpipeline.bam import sortbam, indexbam

import log
logger = log.setup_logger('tagreads',log.get_config())

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
    re.compile( '[A-Z0-9]{5}:\d{1,}:\d{1,}' ),
    re.compile( 'M[0-9]{5}:\d+:\d{9}-[A-Z0-9]{5}:\d:\d{4}:\d{4,5}:\d{4,5}' ),
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
    logger.info( "Gathering existing header for {}".format(bam) )
    hdr = get_rg_headers( bam, SM, CN )
    tag_reads( bam, hdr )

def tag_read( untagged_read, tags ):
    '''
        Tags a given samtools.SamRow with various tags
        Does not replace existing tags

        @param untagged_read - samtools.SamRow
        @param tags - List of valid samspec optional field strings(aka ['RG:Z:Value'])

        @returns a tagged read
    '''
    if untagged_read.FLAG >= 2048:
        # Skip supplementary
        logger.debug( "Skipping read {} because it is supplementary".format(untagged_read.QNAME) )
        return untagged_read
    # Append the new tags
    if untagged_read._tags and untagged_read._tags[-1] != '\t':
        untagged_read._tags += '\t'
    # Append the tags that do not already exist
    # tag+'\t' is used to make it an exact match
    tags = [tag for tag in tags if tag+'\t' not in untagged_read._tags] 
    untagged_read._tags += '\t'.join( tags )
    # Remove trailing tab if all the tags were duplicate
    untagged_read._tags = untagged_read._tags.rstrip()
    # Return the tagged read
    return untagged_read

def tag_readgroup( read ):
    '''
        Tags a given read with the readgroup that the qname
        belongs to depending on how it maps to ID_MAP

        @param read - samtools.SamRow
        
        @returns SamRow that is tagged with the appropriate read group
    '''
    rg = get_rg_for_read( read )
    #logger.debug( "Tagging {} with Read group {}".format(read.qname,rg) )
    return tag_read( read, ['RG:Z:'+rg] )

def tag_reads( bam, hdr ):
    '''
        Sets header of bam and tags all reads appropriately for each platform
        Overwrites existing header
        
        @param bam - Bam file to tag reads in
        @param hdr - Header string to set in the bam(needs newline at the end)
    '''
    # Open the existing bam to fetch the reads to modify from
    untagged_bam = samtools.view( bam )
    # Open a file to write the sam output to with the tagged reads
    samf = bam.replace('.bam','.sam')
    with open( samf, 'w' ) as sam:
        # Write the hdr to the file first
        sam.write( hdr )
        # Tag the reads
        logger.info( "Tagging reads for {}".format(bam) )
        for read in untagged_bam:
            samrow = samtools.SamRow(read)
            read = tag_readgroup( samrow )
            sam.write( str(read) + '\n' )
    # Close stdout
    untagged_bam.close()
    logger.info( "Finished tagging reads for {}".format(bam) )
    logger.info( "Sorting {}".format(bam) )
    b = samtools.view( samf, h=True, S=True, b=True )
    sortbam( b, bam )
    # Close the fh
    b.close()
    # Remove temp sam file
    # maybe some day could even just use pipes all the way through :)
    os.unlink( samf )
    logger.info( "Indexing {}".format(bam) )
    indexbam( bam )

def get_rg_for_read( aread ):
    ''' Gets the read group name for the given samtools.SamRow '''
    rname = aread.QNAME
    for i, p in enumerate( ID_MAP ):
        if p.match( rname ):
            return IDS[i]
    raise UnknownReadNameFormat( "{} is from an unknown platform and cannot be tagged".format(rname) )

def get_rg_headers( bam, SM=None, CN=None ):
    old_header = get_bam_header( bam ) + '\n'

    for id, pl in zip( IDS, PLATFORMS ):
        # Skip headers that exist already
        if 'ID:{}\t'.format(id) in old_header:
            continue

        rg = '@RG\tID:{}\tSM:{}\t'
        if SM is None:
            SM = os.path.basename(bam).replace( '.bam', '' )

        if CN is not None:
            rg += 'CN:{}\tPL:{}'
            old_header += rg.format( id, SM, CN, pl ) + '\n'
        else:
            rg += 'PL:{}'
            old_header += rg.format( id, SM, pl ) + '\n'

    return old_header

def get_bam_header( bam ):
    r = samtools.view( bam, H=True )
    hdr = r.read()
    r.close()
    return hdr.rstrip()

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

