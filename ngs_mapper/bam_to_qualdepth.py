#!/usr/bin/env python2

import sys
import argparse
import json

import bam
import bqd

def main():
    args = parse_args()
    print_json( args )

def print_json( args ):
    # TODO
    # Need to fix this so it uses samtools.mpilup
    pileup = bqd.mpileup( args.bamfile )
    pileup = bqd.parse_pileup( pileup )
    set_unmapped_mapped_reads( args.bamfile, pileup )
    print json.dumps( pileup )

def set_unmapped_mapped_reads( bamfile, pileup ):
    ''' add mapped/unmapped reads to json for each reference '''
    idxstats = bam.get_refstats( bamfile )
    if '*' not in idxstats:
        pileup['unmapped_reads'] = 0
    else:
        pileup['unmapped_reads'] = int(idxstats['*'][3])
    for ref, stats in idxstats.iteritems():
        if ref != '*':
            try:
                pileup[ref]['mapped_reads'] = int(stats[2])
            except KeyError:
                # Skips references that are missing due to no depth
                continue

def parse_refs( bamfile ):
    '''
        Parse all the pileup information for each reference in the BAM file

        @param bamfile - Path to bam file
        
        @returns {refname:{length:,maxdepth:,mindepth:,maxq:,minq:,depth:[],quals:[]}}
    '''
    refs = {}
    for refinfo, ainfo in bam.alignment_info(bamfile):
        maxd,mind,maxq,minq,depths,quals = parse_ainfo( ainfo )
        rname,len,nreads,nureads = refinfo
        refs[rname] = {
            'length':len,
            'maxd':maxd,
            'mind':mind,
            'maxq':maxq,
            'minq':minq,
            'depths':depths,
            'quals':quals
        }

    return refs

def parse_ainfo( ainfo ):
    '''
        Parse the information for each reference
    
        @param ainfo - bam.alignment_info()[1] list

        @returns (maxdepth,mindepth,maxquality,minquality,[depth at each position],[avg quality at each position])
    '''
    maxq = maxd = 0
    minq = mind = sys.maxint
    
    depths = []
    quals = []
    for info in ainfo:
        d,mq,miq,aq = parse_info( info )
        maxq = max( maxq, mq )
        minq = min( minq, miq )
        maxd = max( maxd, d )
        mind = min( mind, d )
        depths.append( d )
        quals.append( aq )

    return (maxd,mind,maxq,minq,depths,quals)

def parse_info( info ):
    '''
        Parse each positions information

        @param info - bam.parse_pileup() output

        @returns (depth,maxquality,minquality,avgquality)
    '''
    pos,b,depth,seq,quals = info
    maxq = 0
    minq = 100
    qsum = 0.0
    for q in quals:
        maxq = max( maxq, q )
        minq = min( minq, q )
        qsum += q
    return (depth,maxq,minq,qsum/len(quals))

def parse_args():
    parser = argparse.ArgumentParser(
        description='Create Quality/Depth json output'
    )

    parser.add_argument(
        dest='bamfile',
        help='Bam file to get gaps for'
    )

    return parser.parse_args()

if __name__ == '__main__':
    main()
