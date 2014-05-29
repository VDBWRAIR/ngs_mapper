#!/usr/bin/env python

import sys
import os
import pprint

# miseqpipeline should be first in the PATH
sys.path.insert( 0, os.environ['PATH'].split( ':' )[0] )
import samtools

def main( args ):
    for col in cols( args.bamfile, args.regionstr, args.minmq, args.minbq ):
        stats = col.base_stats()
        print "{}:".format(col.pos)
        pprint.pprint(stats)

def cols( bamfile, regionstr, *args, **kwargs ):
    ''' Generator for mpileup columns '''
    piles = samtools.mpileup( bamfile, regionstr, args[0], args[1] )
    for pile in piles:
        yield samtools.MPileupColumn( pile )

def parse_args( args=sys.argv[1:] ):
    import argparse
    parser = argparse.ArgumentParser(
        description='''Output mpileup bases and qualities in human readable format'''
    )

    parser.add_argument(
        'bamfile',
        help='Bam file location'
    )

    parser.add_argument(
        'regionstr',
        help='Region string'
    )

    parser.add_argument(
        '-minmq',
        dest='minmq',
        default=20,
        type=int,
        help='Minimum mapping quality to exclude'
    )

    parser.add_argument(
        '-minbq',
        dest='minbq',
        default=25,
        type=int,
        help='Minimum base quality to exclude'
    )

    return parser.parse_args( args )

if __name__ == '__main__':
    main( parse_args() )
