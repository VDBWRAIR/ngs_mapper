#!/usr/bin/env python

## Untested

import argparse
import sys
import itertools
from collections import OrderedDict
import samtools

def main(args):
    return stats_at_pos( args.bamfile, args.regionstr, args.minmq, args.minbq, args.maxd )

def parse_args(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        description='''Gives stats about a given site in a bam file''',
        epilog='''You might use this command to get a list of available reference 
names to use for the regionstr. In the future there will be a list command for this
but for now use this:                                                 
samtools idxstats <in.bam> | awk '!/^*/ {print $1}' | sort | uniq'''
    )

    parser.add_argument(
        dest='bamfile',
        help='Bam file path for stats'
    )

    parser.add_argument(
        dest='regionstr',
        help='Region string to use in format of {refname}:{start}-{stop} where ' \
            'refname should be a reference identifier inside of the bamfile, ' \
            'start and stop are the reference position to look at. Example: ' \
            'Den1Reference:1046-1046. You can see more about regionstr in the ' \
            ' samtools documentation.'
    )

    default_minmq = 25.0
    parser.add_argument(
        '-mmq',
        '--min-mapping-qual',
        dest='minmq',
        type=float,
        default=default_minmq,
        help='Minimum mapping quality to be included in stats. Keep reads that are >= [Default: {}]'.format(default_minmq)
    )

    default_minbq = 25.0
    parser.add_argument(
        '-mbq',
        '--min-base-qual',
        dest='minbq',
        type=float,
        default=default_minbq,
        help='Minimum base quality to be included in stats. Keep bases that are >= [Default: {}]'.format(default_minbq)
    )

    default_maxdepth = 100000
    parser.add_argument(
        '-m',
        '--max-depth',
        dest='maxd',
        type=int,
        default=default_maxdepth,
        help='Maximum read depth at position to use[Default: {}]'.format(default_maxdepth)
    )
    
    return parser.parse_args(args)

def stats_at_pos( bamfile, regionstr, minmq, minbq, maxd ):
    base_stats = compile_stats( stats( bamfile, regionstr, minmq, minbq, maxd ) )
    print "Maximum Depth: {}".format(maxd)
    print "Minumum Mapping Quality Threshold: {}".format(minmq)
    print "Minumum Base Quality Threshold: {}".format(minbq)
    print "Average Mapping Quality: {}".format(base_stats['AvgMapQ'])
    print "Average Base Quality: {}".format(base_stats['AvgBaseQ'])
    print "Depth: {}".format(base_stats['TotalDepth'])
    for base, bstats in base_stats['Bases'].iteritems():
        print "Base: {}".format(base)
        print "\tDepth: {}".format( bstats['Depth'] )
        print "\tAverage Mapping Quality: {}".format( bstats['AvgMapQ'] )
        print "\tAverage Base Quality: {}".format( bstats['AvgBaseQ'] )
        print "\t% of Total: {}".format( bstats['PctTotal'] )

    return base_stats

def stats( bamfile, regionstr, minmq, minbq, maxd ):
    out = samtools.mpileup( bamfile, regionstr, minmq, minbq, maxd )
    
    try:
        o = out.next()
        col = samtools.MPileupColumn( o )
        out.close()
        return col.base_stats()
    except StopIteration:
        return {
            'depth': 0,
            'mqualsum': 0.0,
            'bqualsum': 0.0
        }

def compile_stats( stats ):
    '''
        @param stats - {'depth': 0, 'mqualsum': 0, 'bqualsum': 0, 'ATGCN*..': [quals]} depth is total depth at a position and qualsum is sum of all quality scores bqualsum is read quality sums ATGCN* will be keys for each base seen and the list of quality scores for them
        @return - Dictionary of stats at each base and overall stats {'Bases': {'A': [quals], 'depth': 0, 'avgqual': 0.0}}
    '''
    if stats['depth'] == 0:
        return {
            'TotalDepth': 0,
            'AvgMapQ': 0,
            'AvgBaseQ': 0,
            'Bases': {}
        }
    base_stats = {}
    base_stats['TotalDepth'] = stats['depth']
    base_stats['AvgMapQ'] = round(stats['mqualsum']/stats['depth'],2)
    base_stats['AvgBaseQ'] = round(stats['bqualsum']/stats['depth'],2)
    base_stats['Bases'] = {}
    for base, quals in stats.iteritems():
        # Only interested in base stats in this loop
        if base not in ('depth','mqualsum','bqualsum'):
            if base not in base_stats['Bases']:
                base_stats['Bases'][base] = {}
            mquals = quals['mapq']
            bquals = quals['baseq']
            base_stats['Bases'][base]['Depth'] = len(mquals)
            base_stats['Bases'][base]['AvgMapQ'] = round(float(sum(mquals))/len(mquals),2)
            base_stats['Bases'][base]['AvgBaseQ'] = round(float(sum(bquals))/len(bquals),2)
            base_stats['Bases'][base]['PctTotal'] = round((float(len(mquals))/stats['depth'])*100,2)

    # Quit out of loop we are done
    # Order bases by PctTotal descending
    sorted_bases = sorted( base_stats['Bases'].items(), key=lambda x: x[1]['PctTotal'], reverse=True )
    base_stats['Bases'] = OrderedDict(sorted_bases)
    return base_stats

if __name__ == '__main__':
    main(parse_args())
