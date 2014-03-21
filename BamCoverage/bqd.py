#!/usr/bin/env python

from subprocess import Popen, PIPE
import sys
import json

def main():
    if not len(sys.argv) == 2:
        print "Need bam file"
        sys.exit(-1)

    print json.dumps( parse_pileup( mpileup( sys.argv[1] ) ) )

def mpileup( bamfile, regionstr=None ):
    cmd = ['samtools','mpileup', '-d', '100000']
    if regionstr:
        cmd += ['-r',regionstr]
    cmd.append( bamfile )
    p = Popen( cmd, stdout=PIPE )
    return p.stdout

def parse_pileup( pileup ):
    '''
        Parses the raw pileup output from samtools mpileup and returns a dictionary
        with stats for every reference in the pileup
            maxd/mind - max/min depth found for that reference
            maxq/minq - max/min quality found for that reference
            depths - depth at each base position
            avgquals - average quality at each base position
            length - length of reference

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

if __name__ == '__main__':
    main()
