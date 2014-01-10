#!/usr/bin/env python

import sys
import os
from os.path import *
from subprocess import check_output, PIPE
import argparse
from BamCoverage import bqd, mkg
from BamCoverage.bam_to_qualdepth import set_unmapped_mapped_reads
import json

def main( args ):
    args = handle_args( args )
    jfile = make_json( args.bamfile, args.outpath )
    pngfile = make_image( jfile, args.outpath )

def make_json( bamfile, outpathprefix ):
    pileup = bqd.mpileup( bamfile )
    stats = bqd.parse_pileup( pileup )
    set_unmapped_mapped_reads( bamfile, stats )
    outfile = outpathprefix + '.qualdepth.json'
    with open( outfile, 'w' ) as fh:
        json.dump( stats, fh )

    return outfile

def make_image( jfile, outpathprefix ):
    outfile = outpathprefix + '.qualdepth.png'
    mkg.make_graphic( jfile, outfile, titleprefix=basename(outpathprefix) )
    return outfile
    
def handle_args( args ):
    if args.outprefix is not None:
        outprefix = args.outprefix
    else:
        bamname = basename( args.bamfile )
        outprefix = bamname

    args.outpath = join( args.outdir, outprefix )

    return args

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Runs bam_to_qualdepth.py as well as mkg.py from BamCoverage'
    )

    parser.add_argument(
        'bamfile',
        help='Path to bamfile'
    )

    parser.add_argument(
        '-od',
        default=os.getcwd(),
        dest='outdir',
        help='Where to place the output files(outdir/outprefix.{png,json})[Default: Current directory]'
    )

    parser.add_argument(
        '-op',
        '--out-prefix',
        dest='outprefix',
        default=None,
        help='How to name the output files. Default is to prefix with the bamfile name'
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main(parse_args())
