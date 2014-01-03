#!/usr/bin/env python


## Untested

# the lines to generate the filtered vcf 

# samtools mpileup -uD -f Thai_H.fasta 00006.3_06__RL23__Den1.bam | bcftools view -bvcg - > RAL_samtools.raw.bcf



# bcftools view RAL_samtools.raw.bcf | vcfutils.pl varFilter -D100 > RAL_samtools.vcf


import re
import sys
import subprocess
import os



import argparse
parser = argparse.ArgumentParser(description="Variant caller")
parser.add_argument(
	dest='reference',
	help='Path to reference file'
)
parser.add_argument(
	dest='bam',
	help='Bam file path'
)
'''
parser.add_argument(
	'-d',
	'--depth',
	dest='depth',
	help='Min depth for variant calling'
)
'''

parser.add_argument(
	'-o',
	dest='output',
	help='Output file name',
	default='varient_called_file.vcf'
)
parser.add_argument(
	'--use-filter',
	dest='filteron',
	default=False,
	action='store_true',
	help='Only output variants(use -v option in bcftools)'
)
parser.add_argument(
	'-Q',
	dest='minreadqual',
	default=13,
	help='Minimum read quality to be considered'
)

args = parser.parse_args()


"""
# system error if there are too many arguments entered on the command line
if len(sys.argv) != 3:
    print('Error, please check the number of arguments')
    sys.exit()
"""
# Input the required information

# Reference file
reference_fasta = args.reference 

# Alignment file
alignment_bam = args.bam 

# Minimum Read Quality
minreadqual = args.minreadqual

# Set the threadhold
# threashold = int(sys.argv[3])


# The out put file - must end in .vcf
output_VCF_File = args.output

filter_option = ''
if args.filteron:
	filter_option = 'v'

# call runs an external program and waits for it to quit

from subprocess import call 

for command in ["samtools mpileup -uD -Q {} -f {} {} | bcftools view -cg{} - > {}".format(minreadqual,reference_fasta,alignment_bam,filter_option,output_VCF_File)]:
               #  "bcftools view RAL_samtools.raw.bcf | vcfutils.pl varFilter -D100 > output_VCF_File"):

    # shell=True is so you can handle redirects like in the 3rd command
    print command
    call(command, shell=True)




