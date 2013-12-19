import tempfile
import shutil
import os
from os.path import *

tdir = tempfile.mkdtemp(prefix='tests',suffix='class')

def setUpPackage():
    os.chdir(tdir)

def tearDownPackage():
    shutil.rmtree(tdir)
