#!/usr/bin/env python

from subprocess import Popen, PIPE
import sys
import json

from collections import namedtuple
from itertools import izip

from matplotlib.lines import Line2D

import log
logger = log.setup_logger(__name__, log.get_config())

# Alias our region strings
G = 'Gap'
N = 'Normal'
LC = 'LowCoverage'
LQ = 'LowQuality'
LCQ = 'LowCovQual'

# As a list
REGIONTYPES = [
    G, N, LC, LQ, LCQ
]

def mpileup( bamfile, regionstr=None ):
    '''
    .. todo::
       This needs to go as samtools.mpileup exists
    '''
    cmd = ['samtools','mpileup', '-d', '100000']
    if regionstr:
        cmd += ['-r',regionstr]
    cmd.append( bamfile )
    p = Popen( cmd, stdout=PIPE )
    return p.stdout

def parse_pileup( pileup ):
    '''
    .. todo::
        This needs to be replaced by samtools stuff        

    Parses the raw pileup output from samtools mpileup and returns a dictionary
    with stats for every reference in the pileup

        - maxd/mind - max/min depth found for that reference
        - maxq/minq - max/min quality found for that reference
        - depths - depth at each base position
        - avgquals - average quality at each base position
        - length - length of assembly

    @pileup - file like object that returns lines from samtools mpileup

    @returns dictionary {'ref1': {maxd:0,mind:0,maxq:0,minq:0,depths:[],avgquals:[],length:0}, 'ref2':...}
    '''
    refs = {}
    lastpos = {}
    for line in pileup:
        info = line.rstrip().split('\t')
        # Depth is 0
        if len( info ) == 4:
            ref, pos, n, depth = info
            seq = ''
            quals = ''
        elif len( info ) == 6:
            ref,pos,n,depth,seq,quals = info
        else:
            raise ValueError( "mpileup line {} is unparseable".format(line) )
        # Initialize new reference
        if ref not in refs:
            lastpos[ref] = 0
            refs[ref] = {
                'maxd': 0,
                'mind': 8000,
                'maxq': 0,
                'minq': 99,
                'depths': [],
                'avgquals': [],
                'length': 0
            }
        pos = int(pos)
        # Fill in gaps with blanks
        if lastpos[ref] != pos-1:
            refs[ref]['mind'] = 0
            refs[ref]['minq'] = 0.0
            # From the position 1 past the last
            #  all the way up to the current
            for i in range( lastpos[ref]+1, pos ):
                refs[ref]['depths'].append( 0 )
                refs[ref]['avgquals'].append( 0.0 )
                refs[ref]['length'] += 1
        depth = int(depth)
        refs[ref]['maxd'] = max( refs[ref]['maxd'], depth )
        refs[ref]['mind'] = min( refs[ref]['mind'], depth )
        sum = 0
        for o in quals:
            q = float( ord(o)-33 )
            refs[ref]['maxq'] = max( refs[ref]['maxq'], q )
            refs[ref]['minq'] = min( refs[ref]['minq'], q )
            sum += q
        if depth != 0:
            refs[ref]['avgquals'].append( sum / depth )
        else:
            refs[ref]['avgquals'].append( 0 )
        refs[ref]['depths'].append( depth )
        refs[ref]['length'] += 1
        lastpos[ref] = pos

    return refs

# Named tuple to store each region in
CoverageRegion = namedtuple('CoverageRegion', ['start','end','type'])

def get_region_type(depth, qual, gap, lowqual, lowcov):
    '''
    Return the region type for the given depth and quality combination
    '''
    # Check for gap
    if depth <= gap:
        return G
    # Check for lowcovqual
    elif depth < lowcov and qual < lowqual:
        return  LCQ
    # Check for lowcov
    elif depth < lowcov:
        return LC
    elif qual < lowqual:
        return LQ
    else:
        return N

def regions_from_qualdepth(qualdepth, gap, lowqual, lowcov):
    '''
    Turn qualdepth into a generator of regions
    Each region returned will be a namedtuple with start, end, type filled out
    
    Types are of the following:
      'Gap'
      'Normal'
      'LowCoverage'
      'LowQuality'
      'LowCovQual'
    
    qualdepth is a dictionary with the following keys
        depths(list)
        avgquals(list)
        minq(float)
        maxq(float)
        mind(int)
        maxd(int)
        length(int)
        mapped_reads(int)
        
    lowqual and lowcov define the minimum requirements to be called that
    type of region. This is a non-inclusive operation(aka value < LowCoverage would be called LowCoverage)
    Gap is very similar, except it is inclusive since it could be 0

    Since a visible line requires a start and end point to be separated by 1 unit, each base position is 
    represented via start: baseposition end: baseposition+1
    So a region of

    (1, 2, 'Gap') would represent a gap for only the very first base position
    Where on a reference of length 10;
    (1,11,'Normal') would represent bases 1-10 all normal coverage
    '''
    # The current region we are working on
    curregion = [0,0,'']
    # Loop through depth and avgqualities together
    for basepos, da in enumerate(izip(qualdepth['depths'], qualdepth['avgquals']),start=1):
        # Split up the tuple
        depth, avgqual = da
        # Will hold current region type
        regtype = get_region_type(depth, avgqual, gap, lowqual, lowcov)
        # For the very first baseposition
        if basepos == 1:
            curregion = [basepos,basepos,regtype]
        # If the region type has changed and not the first base
        #  then yield the current region
        elif regtype != curregion[2]:
            # End the current region
            curregion[1] = basepos
            # Yield the built named tuple
            yield CoverageRegion._make(curregion)
            # Start a new region now
            curregion = [basepos,basepos,regtype]
        else:
            curregion[1] = basepos
    # Increment end by 1 to include the end
    curregion[1] += 1
    # reference length +1 because lines require start and end points so
    # bases are start:basepos end:basepos + 1
    if 'reflen' in qualdepth:
        reflen = qualdepth['reflen'] + 1
    else:
        reflen = qualdepth['length'] + 1
    # extend current region if already a gap
    print "Reflen: {0}".format(reflen)
    print "Curregion: {0}".format(curregion)
    if curregion[1] < reflen and curregion[2] == G:
        # Extend as the last region was already a gap
        print 'Extending {0} to end of reference({1})'.format(curregion, reflen)
        curregion[1] = reflen
    # Yield the last build region
    yield CoverageRegion._make(curregion)
    # Build a new gap region at end if needed
    if curregion[1] < reflen:
        region = [curregion[1],reflen,G]
        print 'Creating new gap region to cover end of reference: {0}'.format(region)
        yield CoverageRegion._make(region)

def lines2d_from_regions(yval, regions, line2dargs):
    '''
    Create Line2D's for each region using the start and stop values
    for the x1,x2 values and yval for y1 and y2
    
    line2dargs is a dictionary mapping each region.type to a line2d argument that
    is just all arguments for Line2D to set
    Color the line based on the linecolors provided which are a dictionary
    mapping the the CoverageRegion.type to a color
    
    Returns a generator of Line2D objects
    '''
    for region in regions:
        yield Line2D([region.start,region.end], [yval,yval], **line2dargs[region.type])
