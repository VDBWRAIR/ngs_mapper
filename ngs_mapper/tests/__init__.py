import tempfile
import shutil
import os
from os.path import *
from glob import glob

tdir = tempfile.mkdtemp(prefix='ngs_mappertests',suffix='package')
before = []

def setUpPackage():
    global before
    os.chdir(tdir)
    vpath = join(tdir, '.venv')
    before = glob( '/dev/shm/*' ) + glob( '/tmp/*' )

def tearDownPackage():
    shutil.rmtree(tdir)
    shmcontents = glob( '/dev/shm/*' )
    tmpcontents = glob( '/tmp/*' )

    # Clean up all shm and tmp stuff that was made
    for df in shmcontents + tmpcontents:
        if df not in before:
            print "Removing {}".format(df)
            if isdir( df ):
                shutil.rmtree(df)
            else:
                os.unlink(df)

def setUp( self ):
    super(Base,self).setUp()
