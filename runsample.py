#!/usr/bin/env python

import argparse
import subprocess
import sys
import os

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Runs a single sample through the pipeline'
    )

    parser.add_argument(
        dest='readsdir',
        help='Directory that contains reads to be mapped'
    )

    parser.add_argument(
        dest='reference',
        help='The path to the reference to map to'
    )

    parser.add_argument(
        dest='prefix',
        help='The prefix to put before every output file generated'
    )

    default_outdir = os.getcwd()
    parser.add_argument(
        '-od',
        '--outdir',
        dest='outdir',
        default=default_outdir,
        help='The output directory for all files to be put[Default: {}]'.format(default_outdir)
    )

    return parser.parse_args( args )

# Everything to do with running a single sample
# Geared towards running in a Grid like universe(HTCondor...)

# All the initial Inputs
#READDIR=$1
#REFERENCE=$2

# The idea is that you should only have to run each script almost like just writting in a text file 
#  run_bwa_on_samplename.py READDIR REFERENCE -o BAMFILE
#  varcaller.py BAMFILE REFERENCE -o VARIANTS
#  graphsample.py BAMFILE
