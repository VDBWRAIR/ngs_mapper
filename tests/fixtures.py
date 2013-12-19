import os
from os.path import *

# This directory
THIS=dirname(abspath(__file__))

def ungiz( filepath, dest=os.getcwd() ):
    ''' unpack filepath into dest '''
    import gzip
    fofile = join(dest,basename(filepath).replace('.gz',''))
    with gzip.open(filepath) as fh:
        with open(fofile,'w') as fo:
            fo.write( fh.read() )

    return fofile

def get_sample_paired_reads_compressed():
    path = join(THIS,'fixtures','paired')
    return (join(path,'F.fq.gz'),join(path,'R.fq.gz'),join(path,'ref.fna.gz'))

def get_sample_paired_reads(dir=None):
    paired = get_sample_paired_reads_compressed()
    if not dir:
        import tempfile
        d = tempfile.mkdtemp(prefix='tests')
    else:
        d = dir
    a = ungiz( paired[0], d )
    b = ungiz( paired[1], d )
    c = ungiz( paired[2], d )
    return (a,b,c)
