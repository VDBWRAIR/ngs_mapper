#!/usr/bin/env python2.7

import sys
import os
import subprocess
import re
import os.path
import argparse
import json
import itertools
from StringIO import StringIO

import ngs_mapper.bam
#from bam import parse_pileup, alignment_info, get_refstats
from ngs_mapper.bam import get_refstats
from bam import alignment_info

class RegionTypes(object):
    # It is ok to be lazy sometimes
    n = 'Normal'
    g = 'Gap'
    lc = 'LowCoverage'
    lq = 'LowQuality'

class Thresholds(object):
    # Low Coverage Threshold
    #  Depth at a given alignment position
    #  <= this number is 
    # Low Coverage
    LCDT = 9
    # Gap
    GDT = 0

    # Average Quality Threshold
    #  Average base quality at a given alignment position
    #   <= this number is
    # Low Quality
    LCAQT = 25
    # Gap
    GAQT = 0

    def __str__( self ):
        import pprint
        pprint.pprint( self.__dict__ )

TB = '\t'
SP = ' '
NL = '\n'
INDENT = 16

def main():
    global TB, SP, NL, INDENT
    args = parse_args()

    if args.smalljson:
        TB = SP = NL = ''
        INDENT = None

    comma = ''
    if args.closejson:
        sys.stdout.write( '['+NL )

    for bamfile, samplename in itertools.izip_longest( args.bamfile, args.samplename, fillvalue='' ):
        if not samplename:
            samplename = os.path.basename( bamfile )
        sys.stdout.write( comma )
        output_json( alignment_info( bamfile, args.regionstr ), samplename )
        comma = ','+NL

    if args.closejson:
        sys.stdout.write( ']' )

def output_json( alignment_info, samplename ):
    '''
        Formats output from alignment_info into json
    '''
    sys.stdout.write( TB+'{'+NL )
    sys.stdout.write( TB*2+'"name": "{}",'.format( samplename )+NL ) 
    sys.stdout.write( TB*2+'"references": {'+NL )
    comma = ''
    for refinfo, ai in alignment_info:
        stats = parse_mapstats( ai )
        # Is there a gap at the end that needs to be added
        if stats[-1][1] != refinfo[1]:
            # Is the last region already a gap that we can
            #  add to?
            if stats[-1][2] == RegionTypes.g:
                stats[-1][1] = int(refinfo[1])
            else:
                stats.append( (stats[-1][1], int(refinfo[1]), RegionTypes.g) )
        jsonstr = TB*3+'"' + refinfo[0] + '":{'+NL
        jsonstr += TB*4+'"reflen":' + refinfo[1] + ','+NL
        jsonstr += TB*4+'"mappedreads":' + refinfo[2] + ','+NL
        jsonstr += TB*4+'"regions":' + json.dumps( stats, indent=INDENT ).replace(' ', SP) + NL
        jsonstr += TB*3+'}'
        sys.stdout.write( comma + jsonstr )
        comma = ','+NL
    sys.stdout.write( NL+TB*2+'}'+NL )
    sys.stdout.write( TB+'}' )

def parse_mapstats( mapstats ):
    '''
        Consolidates mapstats from parse_pileup into regions of
            low coverage, gaps and regular

        Beginning of reference is position 1(aka 1 indexed sequences)
        Each adjacent region will share a base. That is, the left region's end == right region's start

        @param mapstats - dictionary of reference information from parse_pileup for single reference
        
        @returns list of regions as tuples (start,stop,regiontype)
    '''
    # Beginning region
    rtype = region_type(mapstats['depths'][0],mapstats['avgquals'][0])
    regions = [(1,0,rtype)]
    for pos,dq in enumerate( zip(mapstats['depths'],mapstats['avgquals']), 1 ):
        depth, avgqual = dq
        rtype = region_type( depth, avgqual )
        # If current region's rtype is different
        # Start a new region
        if regions[-1][2] != rtype:
            # End the last region
            regions[-1] = (regions[-1][0],int(pos),regions[-1][2])
            # Start the new region
            regions.append((pos,0,rtype))

    # Finish up the last region
    #regions[-1][1] = int(pos)
    regions[-1] = (regions[-1][0],int(pos),regions[-1][2])

    return regions

def region_type( depth, avgqual ):
    '''
        Determins the region type from the given depth and average quality
        
        @returns one of the RegionTypes
    '''
    #pos, base, depth, seq, quals = mstat
    #avgqual = sum(quals)/depth
    if depth <= Thresholds.GDT or avgqual <= Thresholds.GAQT:
        return RegionTypes.g
    elif depth <= Thresholds.LCDT:
        return RegionTypes.lc
    elif avgqual <= Thresholds.LCAQT:
        return RegionTypes.lq
    else:
        return RegionTypes.n

def parse_args():
    parser = argparse.ArgumentParser(
        description='Get gap information from indexed BAM file'
    )

    parser.add_argument(
        dest='bamfile',
        nargs='+',
        help='Bam file to get gaps for'
    )

    parser.add_argument(
        '-s',
        '--samplename',
        default=[],
        nargs='+',
        help='Name of sample. Default is to use bamfile name'
    )

    parser.add_argument(
        '--close-json',
        dest='closejson',
        default=False,
        action='store_true',
        help='Should the json generated have an opening [ and closing ]'
    )

    parser.add_argument(
        '--small-json',
        dest='smalljson',
        default=False,
        action='store_true',
        help='Should the json output be compressed(no whitespaces)'
    )

    parser.add_argument(
        '-r',
        '--region',
        dest='regionstr',
        default=None,
        help='samtools region to run on. Default' \
            ' is to run on all chr & regions.' \
            'Format is: "chr":start-stop where' \
            ' chr is maybe a reference name' \
            ' start & stop are 1 indexed positions along chr'
    )

    return parser.parse_args()

if __name__ == '__main__':
    main()
