#!/usr/bin/env python


## Untested

# the lines to generate the filtered vcf 



import re
import sys
import subprocess
import os


# Arguments that need to be entered into the command line
import argparse
parser = argparse.ArgumentParser(description="Variant caller")

parser.add_argument(
	dest='bam',
	help='Bam file path'
)


parser.add_argument(
	'-o',
	dest='output',
	help='Output file name',
	default='varient_called_file.vcf'
)


parser.add_argument(
	dest='reference',
	help='Path to reference file'
)


parser.add_argument(
	'-b',
	dest='Alt_Ref_Ratio',
	default=0.19,
	help='Require the ratio of alt/ref above a specified threashold'
)


parser.add_argument(
	'-bq',
	dest='base_quality',
	default=20,
	help='Require the base quality threashold'
)


parser.add_argument(
	'-mq',
	dest='map_quality',
	default=25,
	help='Require the map quality threashold'
)


parser.add_argument(
	'-s',
	dest='strand_bias',
	default=0.0001,
	help='Require the map strand bias threashold'
)


parser.add_argument(
	'-a',
	dest='reads',
	default=10,
	help='Require the least number of reads supporting each strand for alternative allele'
)




"""
# system error if there are too many arguments entered on the command line
if len(sys.argv) != 3:
    print('Error, please check the number of arguments')
    sys.exit()
"""

args = parser.parse_args()


# Input the required information


# Alignment file
alignment_bam = args.bam 

# The out put file - must end in .vcf
output_VCF_File = args.output

# Reference file
ref_fa = args.reference 

# Alt:Ref file
alt_ref = args.Alt_Ref_Ratio

# Base Quality file
base = args.base_quality 

# Map Quality file
mapq = args.map_quality

# Strand Bias file
strand = args.strand_bias

# Supporting reads file
read = args.reads

'''
filter_option = ''
if args.filteron:
	filter_option = 'v'
'''


# call runs an external program and waits for it to quit

from subprocess import call 

# java jar SNVerIndividual -i xxxx.bam -o output_vcf -r xxxx.fa -b 0.19 -bq 20 -mq 25 -s 0.0001 -a 10

command = "SNVer -i {} -o {} -r {} -b {} -bq {} -mq {} -s {} -a {}".format(alignment_bam,output_VCF_File,ref_fa,alt_ref,base,mapq,strand,read)
# shell=True is so you can handle redirects like in the 3rd command
print command
call(command, shell=True)
