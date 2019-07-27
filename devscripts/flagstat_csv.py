import os
import sys
from os.path import *
import re
from glob import glob

def parse_flagstats( flagstatfile ):
    '''
    535812 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 duplicates
    521761 + 0 mapped (97.38%:-nan%)
    535812 + 0 paired in sequencing
    267930 + 0 read1
    267882 + 0 read2
    516905 + 0 properly paired (96.47%:-nan%)
    521068 + 0 with itself and mate mapped
    693 + 0 singletons (0.13%:-nan%)
    162 + 0 with mate mapped to a different chr
    162 + 0 with mate mapped to a different chr (mapQ>=5)

    >>> r = parse_flagstats( 'Projects/1514A00405183N-F/flagstats.txt' )
    >>> r['qcpassed']
    535812.0
    >>> r['qcfailed']
    0.0
    >>> r['mapped']
    521761.0
    >>> r['mappedpct']
    97.38
    >>> r['paired']
    535812.0
    >>> r['read1']
    267930.0
    >>> r['read2']
    267882.0
    >>> r['properpair']
    516905.0
    >>> r['properpairpct']
    96.47
    >>> r['mated']
    521068.0
    >>> r['singletons']
    693.0
    >>> r['singletonspct']
    0.13
    >>> r['wrongmatechr']
    162.0
    '''
    p = '(\d+) \+ (\d+) (.*)'
    p2 = '\((\d+\.\d+)%'
    with open( flagstatfile ) as fh:
        contents = fh.read()
        matches = re.findall( p, contents )
    from collections import OrderedDict
    stats = OrderedDict(
        qcpassed = matches[0][0],
        qcfailed = matches[0][1],
        mapped = matches[2][0],
        mappedpct = re.search( p2, matches[2][2] ).group(1),
        paired = matches[3][0],
        read1 = matches[4][0],
        read2 = matches[5][0],
        properpair = matches[6][0],
        properpairpct = re.search( p2, matches[6][2] ).group(1),
        mated = matches[7][0],
        singletons = matches[8][0],
        singletonspct = re.search( p2, matches[8][2] ).group(1),
        wrongmatechr = matches[9][0],
    )
    return {k:float(v) for k,v in stats.items()}

def csv_out( filepath, stats ):
    samplename = basename( dirname( filepath ) )
    return '{},{}\n'.format(samplename,','.join([str(v) for v in stats.values()]))

def main( ):
    flagstats = glob( 'Projects/*/flagstats*' )
    from StringIO import StringIO
    output = StringIO()
    headers = []
    for i, fs in enumerate(flagstats):
        stats = parse_flagstats( fs )
        if headers:
            assert headers == stats.keys(), 'Headers are not in the same order {} != {}'.format(headers,stats.keys())
        headers = stats.keys()
        output.write( csv_out( fs, stats ) )
    print 'samplename,{}'.format(','.join(headers))
    print output.getvalue()

if __name__ == '__main__':
    main()
