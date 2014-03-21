import subprocess
import os.path
import os
import glob
import gzip
import tempfile
import shutil

def get_bwa_path( ):
    ''' Run which command to get bwa path '''
    try:
        return subprocess.check_output( ['which', 'bwa'] ).strip()
    except subprocess.CalledProcessError as e:
        return ''

# Some useful stuff
this_dir = os.path.dirname( os.path.abspath( __file__ ) )
INPUT_PATH = os.path.join( this_dir, 'input.fastq.gz' )
INPUT_SFF_PATH = os.path.join( this_dir, 'input.sff.gz' )
REF_PATH = os.path.join( this_dir, 'ref.fa' )
BWA_PATH = get_bwa_path()

def test_bwa_available( ):
    ''' Make sure bwa is in path and available '''
    assert get_bwa_path != ''

def ungzip( filepath, outpath=None ):
    ''' TODO: Should check to make sure it needs to actually ungzip a file '''
    if outpath is None:
        bn, ext = os.path.splitext( filepath )
    else:
        bn = outpath
    with open( bn, 'wb' ) as fo:
        with gzip.open( filepath ) as fi:
            fo.write( fi.read() )
    return bn

def unpack_files( ):
    fastq = ungzip( INPUT_PATH )
    sff = ungzip( INPUT_SFF_PATH )
    ref = REF_PATH
    return (fastq, sff, ref)

    
class Base( object ):
    @classmethod
    def setUpClass( self ):
        self.tempdir = tempfile.mkdtemp()
        os.chdir( self.tempdir )

    @classmethod
    def tearDownClass( self ):
        os.chdir( '/' )
        shutil.rmtree( self.tempdir )


def create_fakefasta( filename, readno ):
    ''' Create readno fasta sequences '''
    with open( filename, 'w' ) as fh:
        for i in range( readno ):
            fh.write( '>seq{}\nATGC\n'.format(i) )
