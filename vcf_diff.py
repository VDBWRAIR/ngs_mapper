#!/usr/bin/env python

import vcf
import sys

with open(sys.argv[1]) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith( '#' ):
            continue
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1090-01.bam
        try:
            chrm,pos,id,ref,alt,qual,filter,info = line.split('\t')
            info = dict( [i.split('=') for i in info.split(';')] )
            cb = info['CB']
            if cb != ref:
                print "{}\t{}\t{}\t{}".format(
                    pos, ref, alt, info['CB']
                )
        except ValueError as e:
            print line.split( '\t' )
            print e
