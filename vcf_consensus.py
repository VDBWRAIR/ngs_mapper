#!/usr/bin/env python

import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import vcf
import os

def main( args ):
    seqs = iter_refs( args.vcffile )
    write_fasta( seqs, args.output_file )

def write_fasta( records, outputfile ):
    '''
        Writes SeqRecord objects to fasta outputfile

        @param records - Bio.SeqRecord list
        @param outputfile - File path to write output to

        @returns number of records written
    '''
    try:
        numwrote = SeqIO.write( records, outputfile, 'fasta' )
    except AttributeError as e:
        if 'NoneType' not in e.message:
            raise e
        numwrote = 0
    # Can't check to see if records is empty since it may be a generator
    if numwrote == 0:
        os.unlink( outputfile )

    return numwrote

def iter_refs( vcffile ):
    '''
        Iterates over a given vcf file and outputs a list of Bio.Seq.Seq objects
        that represent the consensus sequence for each of the references in the vcffile

        @param vcffile - path to a vcf file

        @returns list of Bio.Seq.Seq objects for every reference in vcffile
    '''
    # Last reference seen
    lastref = ''
    # Current consensus sequence
    consensus = ''
    for row in vcf.Reader( open(vcffile) ):
        # First iteration
        if lastref == '':
            lastref = row.CHROM
        # New Ref so yield our seq
        if row.CHROM != lastref:
            yield SeqRecord(
                Seq( consensus, generic_dna, ),
                id=lastref,
                description='',
                name=''
            )
            lastref = row.CHROM
            consensus = ''

        # Add to the consensus
        consensus += row.INFO['CB']

    yield SeqRecord(
        Seq( consensus, generic_dna, ),
        id=lastref,
        description='',
        name=''
    )

def parse_args( args=sys.argv[1:] ):
    import argparse
    parser = argparse.ArgumentParser(
        description='Creates a consensus fasta file from a base_caller.py gnerated ' \
            'vcf file'
    )

    parser.add_argument(
        'vcffile',
        help='VCF file path to generate a consensus from'
    )

    parser.add_argument(
        '-o',
        dest='output_file',
        default=None,
        help='Output file path for fasta file[Default: name of vcffile with .vcf ' \
            'replaced with .fasta'
    )

    pa = parser.parse_args( args )

    if pa.output_file is None:
        pa.output_file = pa.vcffile.replace('.vcf','.fasta')

    return pa


if __name__ == '__main__':
    main( parse_args() )
