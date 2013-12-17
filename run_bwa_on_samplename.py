import os
from os.path import *
from glob import glob
import sys
import argparse

def main( args ):

def parse_args( argv=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='''Runs all found runs for a samplename against a reference using bwa'''
    )

    return parser.parse_args(argv)

if __name__ == '__main__':
    main(parse_args())
