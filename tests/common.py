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

import fixtures
from bam import indexbam
class BaseBamRef(BaseClass):
    bam = join(fixtures.THIS,'fixtures','varcaller','paired.bam.gz')
    ref = join(fixtures.THIS,'fixtures','varcaller','ref.fasta.gz')
    mytempdir = ''

    @classmethod
    def setUpClass(klass):
        # Unpacks everything once so it doesn't slow down so much
        super(BaseBamRef,klass).setUpClass()
        import tempfile
        klass.mytempdir = tempfile.mkdtemp(prefix='basebamref',suffix='test')
        klass.bam = fixtures.ungiz(klass.bam,klass.mytempdir)
        klass.ref = fixtures.ungiz(klass.ref,klass.mytempdir)
        klass.bamindex = indexbam( klass.bam )

    @classmethod
    def tearDownClass(klass):
        super(BaseBamRef,klass).tearDownClass()
        import shutil
        shutil.rmtree(klass.mytempdir)

class BaseBaseCaller(BaseClass):
    def setUp( self ):
        super( BaseBaseCaller, self ).setUp()
        fixpath = join( fixtures.THIS, 'fixtures', 'base_caller' )
        self.bam = join( fixpath, 'test.bam' )
        self.bai = join( fixpath, 'test.bam.bai' )
        self.ref = join( fixpath, 'testref.fasta' )
        self.sam = join( fixpath, 'test.sam' )
        self.vcf = join( fixpath, 'test.vcf' )
        self.template = join( fixpath, 'template.vcf' )
