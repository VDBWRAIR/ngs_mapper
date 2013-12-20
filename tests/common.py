import tempfile
import shutil
import os
from os.path import *
from . import tdir

class BaseClass( object ):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        self.tempdir = tempfile.mkdtemp(prefix='unit',suffix='test',dir=tdir)
        os.chdir(self.tempdir)

    def tearDown(self):
        os.chdir(tdir)

    def create_file(self,filepath,contents):
        linecount = 0
        with open(filepath,'w') as fh:
            for line in contents.splitlines(True):
                fh.write(line)
                linecount += 1
        return linecount
