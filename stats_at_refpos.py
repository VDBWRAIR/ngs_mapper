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

    parser.add_argument(
        '-mmq',
        '--min-mapping-qual',
        dest='minmq',
        default=25,
        help='Minimum mapping quality to be included in stats. That is, keep reads that are >= [Default: 25]'
    )

    parser.add_argument(
        '-mbq',
        '--min-base-qual',
        dest='minbq',
        default=20,
        help='Minimum base quality to be included in stats[Default: 20]'
    )

    parser.add_argument(
        '-m',
        '--max-depth',
        dest='maxd',
        default=100000,
        help='Maximum read depth at position to use[Default: 100000]'
    )
    
    return parser.parse_args(args)

def mpileup_pysam( bamfile, regionstr, minmq=25, maxd=10000 ):
    import pysam
    samfile = pysam.Samfile( bamfile )
    #samfile.pileup( region=regionstr )
    reference,reg = regionstr.split(':')
    start,stop = reg.split('-')
    rng = range(int(start),int(stop)+1)
    for col in samfile.pileup(region=regionstr):
         if col.pos in rng:
            yield pysam_col( col.pileups, reference, minmq )

def pysam_col( pysamcol, reference, minmq ):
    '''
        @param pysamcol - pysam.PileupProxy.pileups
        @param reference - Reference(chrm) name
        @param minmq - Filter reads with minmq < minmq
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
        print mapq
        # Ignore any reads with minmq below minmq
        if mapq < minmq:
            print 'skipping'
            continue
        # Sometimes mapq is a character and sometimes an int??
        if type(mapq) in (int,long):
            mapq = chr(mapq)
        mapquals.write(mapq)
        readquals.write(aread.qual[qpos])
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
    return mpileup_pysam( bamfile, regionstr, minmq, maxd )

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
    print "Average Base Quality: {}".format(base_stats['AvgReadQ'])
    print "Depth: {}".format(base_stats['TotalDepth'])
    for base, bstats in base_stats['Bases'].iteritems():
        print "Base: {}".format(base)
        print "\tDepth: {}".format( bstats['Depth'] )
        print "\tAverage Mapping Quality: {}".format( bstats['AvgMapQ'] )
        print "\tAverage Read Quality: {}".format( bstats['AvgReadQ'] )
        print "\t% of Total: {}".format( bstats['PctTotal'] )

    return base_stats

def stats( bamfile, regionstr, minmq, minbq, maxd ):
    base_stats = {}
    pile = mpileup( bamfile, regionstr, minmq=minmq, maxd=maxd )
    for line in pile:
        line = line.rstrip().split()
        pos = int(line[1])
        depth = int(line[3])
        bases = []
        if depth > 0:
            print line[6]
            quals = [ord(q) for q in line[6]]
            for i in range(len(line[4])):
                b = line[4][i]
                if b.upper() not in 'ATGC*N':
                    #print "Skipping unknown base '{}'".format(b)
                    # Do insert/delete
                    if b in '+-':
                        quals = indel( line[4], i, quals )
                    continue
                bases.append(b)
        else:
            quals = []
        bq = [x for x in itertools.izip_longest(bases,quals, fillvalue='!')]
        assert len(quals) == len(bases), "Number of quals {} !=  number of bases {}\n{}".format(len(quals),len(bases),line)

        # Stats
        stats = {
            'depth': 0,
            'mqualsum': 0.0,
            'rqualsum': 0.0
        }
        for base, qual in bq:
            base = base.upper()
            if base not in stats:
                stats[base] = []
            stats[base].append( qual )
            stats['depth'] += 1
            stats['mqualsum'] += float(qual)

        return stats

def compile_stats( stats ):
    '''
        @param stats - {'depth': 0, 'mqualsum': 0, 'rqualsum': 0, 'ATGCN*..': [quals]} depth is total depth at a position and qualsum is sum of all quality scores rqualsum is read quality sums ATGCN* will be keys for each base seen and the list of quality scores for them
        @return - Dictionary of stats at each base and overall stats {'Bases': {'A': [quals], 'depth': 0, 'avgqual': 0.0}}
    '''
    if stats['depth'] == 0:
        return {
            'TotalDepth': 0,
            'AvgMapQ': 0,
            'AvgReadQ': 0,
            'Bases': {}
        }
    base_stats = {}
    base_stats['TotalDepth'] = stats['depth']
    base_stats['AvgMapQ'] = round(stats['mqualsum']/stats['depth'],2)
    base_stats['AvgReadQ'] = round(stats['rqualsum']/stats['depth'],2)
    base_stats['Bases'] = {}
    for base, quals in stats.iteritems():
        if base not in ('depth','mqualsum','rqualsum'):
            if base not in base_stats['Bases']:
                base_stats['Bases'][base] = {}
            base_stats['Bases'][base]['Depth'] = len(quals)
            base_stats['Bases'][base]['AvgMapQ'] = round(float(sum(quals))/len(quals),2)
            base_stats['Bases'][base]['AvgReadQ'] = round(0.0,2)
            base_stats['Bases'][base]['PctTotal'] = round((float(len(quals))/stats['depth'])*100,2)

    # Quit out of loop we are done
    # Order bases by PctTotal descending
    sorted_bases = sorted( base_stats['Bases'].items(), key=lambda x: x[1]['PctTotal'], reverse=True )
    base_stats['Bases'] = OrderedDict(sorted_bases)
    return base_stats

if __name__ == '__main__':
    main(parse_args())
