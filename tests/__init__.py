import tempfile
import shutil
import os
from os.path import *

tdir = tempfile.mkdtemp(prefix='tests',suffix='class')

def setUpPackage():
    import test_install
    os.chdir(tdir)
    vpath = join(tdir, '.venv')

def tearDownPackage():
    shutil.rmtree(tdir)
