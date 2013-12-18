import tempfile
import shutil
import os

class BaseClass( object ):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp(prefix='test')
        os.chdir(self.tempdir)

    def tearDown(self):
        os.chdir('/')
        shutil.rmtree(self.tempdir)

