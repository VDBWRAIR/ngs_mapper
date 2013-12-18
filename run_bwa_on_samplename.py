import os
from os.path import *
from glob import glob
import sys
import argparse

from .data import *

READSBYSAMPLEDIR='/home/EIDRUdata/NGSData/ReadsBySample'

def main( args ):

def parse_args( argv=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='''Runs all found runs for a samplename against a reference using bwa'''
    )

    parser.add_argument(
        dest='samplename',
        help='The samplename to map data for'
    )

    parser.add_argument(
        dest='reference',
        help='The reference to map the data to'
    )

    parser.add_argument(
        '-r',
        '--reads-by-sample',
        dest='readsbysampledir',
        default=READSBYSAMPLEDIR,
        help='What is the base path to the ReadsBySample directory'
    )

    return parser.parse_args(argv)

if __name__ == '__main__':
    main(parse_args())
