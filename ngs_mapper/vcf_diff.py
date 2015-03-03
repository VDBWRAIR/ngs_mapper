import vcf
import sys
import argparse

def main():
    args = parse_args()
    vcf_reader = vcf.Reader(open(args.vcf_file, 'r'))
    print 'Reference\tPosition\tReference Base\tCalled Base'
    for record in vcf_reader:
        cb = record.INFO['CB']
        ref = record.REF
        pos = record.POS
        ref_seq = record.CHROM

        if ref != cb:
            print "{}\t{}\t{}\t{}".format(
                ref_seq, pos, ref, cb 
                )

def parse_args( ):
    parser = argparse.ArgumentParser(description='vcf_diff')
    parser.add_argument(
        dest="vcf_file",
        help="VCF File"
    )
    args = parser.parse_args()
    return args
