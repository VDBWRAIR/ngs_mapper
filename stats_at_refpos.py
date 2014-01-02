#!/usr/bin/env python

## Untested

from subprocess import Popen, PIPE
import argparse
import sys
import itertools

def main(args):
    stats_at_pos( args.refpos, args.bamfile, args.minqual, args.maxd )

def parse_args(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        dest='bamfile',
        help='Bam file path for stats'
    )

    parser.add_argument(
        dest='refpos',
        type=int,
        help='Position on reference to get stats for'
    )

    parser.add_argument(
        '-Q',
        '--min-qual',
        dest='minqual',
        default=25,
        help='Minimum read quality to be included in stats'
    )

    parser.add_argument(
        '-m',
        '--max-depth',
        dest='maxd',
        default=100000,
        help='Maximum read depth at position to use'
    )
    
    return parser.parse_args(args)

def mpileup( bamfile, regionstr=None, minqual=25, maxd=50000 ):
    cmd = ['samtools','mpileup','-Q','{}'.format(minqual),'-d','{}'.format(maxd)]
    if regionstr:
        cmd += ['-r',regionstr]
    cmd.append( bamfile )
    p = Popen( cmd, stdout=PIPE )
    return p.stdout

def stats_at_pos( refpos, bamfile, minqual, maxd ):
    base_stats = {}
    pile = mpileup( bamfile )
    for line in pile:
        line = line.rstrip().split()
        pos = int(line[1])
        depth = int(line[3])
        quals = [ord(q)-33 for q in line[5]]
        bases = [b for b in line[4] if b.upper() in 'ATGC' ]
        if pos != refpos:
            continue
        bq = [x for x in itertools.izip_longest(bases,quals, fillvalue='!')]
        assert len(quals) == len(bases), "{} != {}\n{}".format(len(quals),len(bases),bq)

        # Stats
        stats = {
            'depth': 0,
            'avgqual': 0.0
        }
        for base, qual in bq:
            base = base.upper()
            if base not in stats:
                stats[base] = []
            stats[base].append( qual )
            stats['depth'] += 1
            stats['avgqual'] += float(qual)

        print "Minumum Quality Threshold: {}".format(minqual)
        print "Maximum Depth: {}".format(maxd)
        print "Average Quality: {}".format(stats['avgqual']/stats['depth'])
        for base, quals in stats.iteritems():
            if base not in ('depth','avgqual'):
                print "Base: {}".format(base)
                print "\tDepth: {}".format(len(quals))
                print "\tAverage Quality: {}".format(float(sum(quals)) / len(quals))
                print "\t% of Total: {}".format((float(len(quals))/stats['depth'])*100)

        # Quit out of loop we are done
        break

if __name__ == '__main__':
    main(parse_args())
