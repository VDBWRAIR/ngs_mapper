#!/usr/bin/env python

## Untested

from subprocess import Popen, PIPE
import argparse
import sys
import itertools
from collections import OrderedDict

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

    default_minbq = 20.0
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

def mpileup_pysam( bamfile, regionstr, minmq=25, minbq=20, maxd=10000 ):
    import pysam
    samfile = pysam.Samfile( bamfile )
    #samfile.pileup( region=regionstr )
    reference,reg = regionstr.split(':')
    start,stop = reg.split('-')
    rng = range(int(start),int(stop)+1)
    for col in samfile.pileup(region=regionstr):
         if col.pos in rng:
            yield pysam_col( col.pileups, reference, minmq, minbq )

def pysam_col( pysamcol, reference, minmq, minbq ):
    '''
        @param pysamcol - pysam.PileupProxy.pileups
        @param reference - Reference(chrm) name
        @param minmq - Filter reads with minmq < minmq
        @param minbq - Filter bases with minbq < minbq
        @returns 'chrm pos N depth readbases readquals mapquals'
    '''
    from StringIO import StringIO
    depth = 0
    readseq = StringIO()
    readquals = StringIO()
    mapquals = StringIO()
    # Position on reference for the column we are looking at
    refpos = pysamcol[0].alignment.pos + pysamcol[0].qpos
    for pread in pysamcol:
        qpos = pread.qpos-1
        aread = pread.alignment
        mapq = aread.mapq
        bq = aread.qual[qpos]
        # Ignore any reads with minmq below minmq
        if mapq < minmq:
            #print 'skipping mq'
            continue
        #print "{} < {}?".format(ord(bq)-33,minbq)
        if ord(bq)-33 < minbq:
            #print 'skipping bq'
            continue
        # Sometimes mapq is a character and sometimes an int??
        if type(mapq) in (int,long):
            mapq = chr(mapq)
        mapquals.write(mapq)
        readquals.write(bq)
        seq = aread.seq[qpos]
        readseq.write( aread.seq[qpos] )
        depth += 1

    return '\t'.join([
        reference,
        str(refpos),
        'N',
        str(depth),
        readseq.getvalue(),
        readquals.getvalue(),
        mapquals.getvalue()
    ])

def mpileup_popen( bamfile, regionstr=None, minqual=25, maxd=100000 ):
    cmd = ['samtools','mpileup','-Q','{}'.format(minqual),'-d','{}'.format(maxd)]
    if regionstr:
        cmd += ['-r',regionstr]
    cmd.append( bamfile )
    p = Popen( cmd, stdout=PIPE )
    return p.stdout

def mpileup( bamfile, regionstr=None, minmq=25, minbq=20, maxd=100000 ):
    return mpileup_pysam( bamfile, regionstr, minmq, minbq, maxd )

def indel( sequence, pos, quallist, qualpadding=0 ):
    ''' 
        Returns new quallist with inserted quals
    
        >>> from nose.tools import eq_
        >>> quals = [1,1,1,1]
        >>> seq = 'abc+2aag'
        >>> r = indel( seq, 3, quals )
        >>> e = [1,1,1,0,0,1]
        >>> eq_( e, r, "Insert in middle {} != {}".format(e,r) )
        >>> seq = '+2aaAAAA'
        >>> r = indel( seq, 0, quals )
        >>> e = [0,0,1,1,1,1]
        >>> eq_( e, r, "Insert at beginning {} != {}".format(e,r) )
        >>> seq = 'AAAA+2aa'
        >>> e = [1,1,1,1,0,0]
        >>> r = indel( seq, 4, quals )
        >>> eq_( e, r,"Insert at end {} != {}".format(e,r) )
    '''
    # Should be how many as pos is the index of the +/- in sequence
    n = int(sequence[pos+1])
    if sequence[pos] == '-':
        raise ValueError( "I don't know how to do deletions yet" )
    elif sequence[pos] == '+':
        left = quallist[:pos]
        right = quallist[pos:]
        insert = [qualpadding]*n
        quallist = left + insert + right

    return quallist

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
    base_stats = {}
    pile = mpileup( bamfile, regionstr, minmq=minmq, minbq=minbq, maxd=maxd )
    for line in pile:
        line = line.rstrip().split('\t')
        pos = int(line[1])
        depth = int(line[3])
        bases = []
        mapq = []
        baseq = []
        if depth > 0:
            mapq = [ord(q) for q in line[6]]
            baseq = [ord(q)-33 for q in line[5]]
            bases = [b for b in line[4]]
        bq = [x for x in itertools.izip_longest(bases, mapq, baseq, fillvalue='!')]
        lmq = len(mapq)
        lb = len(bases)
        lbq = len(baseq)
        if lmq != lb and lmq != lbq:
            sys.stderr.write( "Number of mapquals {} basequals {} !=  number of bases {}\n".format(lmq, lbq, lb) )
            sys.stderr.write( "Line that caused error:\n{}\n".format(line) )
            #sys.stderr.write( "(base,mapq,baseq) all together:\n{}\n".format(bq) )
            sys.exit( -1 )

        # Stats
        stats = {
            'depth': 0,
            'mqualsum': 0.0,
            'bqualsum': 0.0
        }
        for base, mq, bq in bq:
            base = base.upper()
            if base not in stats:
                stats[base] = {'baseq':[], 'mapq': []}
            stats[base]['mapq'].append( mq )
            stats[base]['baseq'].append( bq )
            stats['depth'] += 1
            stats['mqualsum'] += float(mq)
            stats['bqualsum'] += float(bq)

        return stats

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
