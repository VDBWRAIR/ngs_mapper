import os
from os.path import *
from . import tdir

# This directory
THIS=dirname(abspath(__file__))
FIXDIR=join(THIS,'fixtures')

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
        d = tempfile.mkdtemp(prefix='tests',dir=tdir)
    else:
        d = dir
    a = ungiz( paired[0], d )
    b = ungiz( paired[1], d )
    c = ungiz( paired[2], d )
    return (a,b,c)

def unpack_integrated( where ):
    nppath = join(THIS,'fixtures','nonpaired')
    ppath = join(THIS,'fixtures','paired')
    ipath = join(THIS,'fixtures','integrated')

    readdir = join(where,'reads')
    os.mkdir(readdir)

    refpath = join(ppath,'ref.fna')
    ref = ungiz( refpath + '.gz', where )

    np = join(nppath,os.listdir(nppath)[0].replace('.gz',''))
    np = ungiz( np + '.gz', readdir )

    p1 = join(ppath,'F.fq')
    p1 = ungiz( p1 + '.gz', readdir )
    os.rename( p1, join(dirname(p1),'sample1_S01_L001_R1_001_1979_01_01.fastq') )
    p2 = join(ppath,'R.fq')
    p2 = ungiz( p2 + '.gz', readdir )
    os.rename( p2, join(dirname(p2),'sample1_S01_L001_R2_001_1979_01_01.fastq') )

    from glob import glob
    expected = {basename(file).replace('.gz',''):ungiz(file,where) for file in glob( join(ipath,'*.gz') )}

    expected['NP'] = np
    expected['P'] = (p1,p2)
    expected['REF'] = ref

    return expected
