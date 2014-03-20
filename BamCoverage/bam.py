import sys
import re
from cStringIO import StringIO
from subprocess import Popen, PIPE

def get_refstats( bamfile ):
    '''
        Run samtools idxstats on the given bamfile

        @returns dictionary keyed by refname and values of [refname, reflen, #mapped reads, singletons]
    '''
    cmd = ['samtools','idxstats',bamfile]
    p = Popen( cmd, stdout=PIPE )
    sout,serr = p.communicate()
    return {line.split()[0]:line.split() for line in sout.splitlines()}

def alignment_info( bam, regionstr=None ):
    '''
        Gets information for a bam file for every column in its alignment
        This is slow. It seems that the samtools mpileup is slow too so I
        guess it is just as slow as that so that is good?

        Deprecated: Try using bqd.parse_pileup instead which should be faster but yields
            a slightly different output

        @param bam - Path to bamfile
        @param regionstr - Region string format specified by
            samtools mpileup -r option

        @yields (get_refstat output for current ref, [parse_pileup() dictionary for each reference])
    '''
    refstats = get_refstats( bam )
    from bqd import mpileup, parse_pileup
    pileup = parse_pileup( mpileup( bam, regionstr ) )
    for ref, stats in pileup.iteritems():
        yield refstats[ref], stats
    
    '''
    start = 0
    end = 9999999
    if regionstr is None:
        cols = samfile.pileup()
    else:
        m = re.match( '(.*?):(\d+)-(\d+)', regionstr )
        if not m:
            sys.stdout.write( "Incorrect regionstr given {}\n".format(regionstr) )
            sys.exit(1)
        chr, start, end = m.groups()
        start = int(start)
        end = int(end)
        cols = samfile.pileup( chr, start, end )

    current_ref = None
    ref_mapstats = []
    for pileup in cols:
        # Only process columns between start and end
        if pileup.pos >= start and pileup.pos < end:
            # Fetch the current ref name
            ref = samfile.getrname( pileup.tid )
            # First iteration only
            if current_ref is None:
                current_ref = ref
            # New ref found so yield prior results
            elif current_ref != ref:
                yield (refstats[current_ref], ref_mapstats)
                ref_mapstats = []
                current_ref = ref
            ref_mapstats.append( parse_pileup( pileup ) )
    yield (refstats[current_ref], ref_mapstats)
    '''


def parse_pileup1( pileup ):
    '''
        Iterate through all pileups and return
            dictionary

        @param pileups - pysam.IteratorColumn

        @returns (pos,base,depth,seq,qual)
    '''
    seq = StringIO()
    quals = []
    for pread in pileup.pileups:
        seq.write( pread.alignment.seq[pread.qpos] )
        quals.append( ord(pread.alignment.qual[pread.qpos])-33 )
    return (pileup.pos, 'N', pileup.n, seq.getvalue(), quals)
