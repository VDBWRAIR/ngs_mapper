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
