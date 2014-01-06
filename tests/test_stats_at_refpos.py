import fixtures
import common
from bam import indexbam

from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import MagicMock, patch, Mock, call

from os.path import *
import os

class Base(common.BaseClass):
    bam = join(fixtures.THIS,'fixtures','integrated','paired.bam.gz')
    mytempdir = ''

    @classmethod
    def setUpClass(klass):
        super(Base,klass).setUpClass()
        import tempfile
        klass.mytempdir = tempfile.mkdtemp()
        klass.bam = fixtures.ungiz(klass.bam,klass.mytempdir)
        klass.bamindex = indexbam( klass.bam )

    @classmethod
    def tearDownClass(klass):
        super(Base,klass).tearDownClass()
        import shutil
        shutil.rmtree(klass.mytempdir)

    def setUp( self ):
        super(Base,self).setUp()
        self.bam = Base.bam

class TestMpileup(Base):
    def _CM( self, bamfile, regionstr, minqual, maxd ):
        from stats_at_refpos import mpileup
        return mpileup( bamfile, regionstr, minqual, maxd )

    @patch('stats_at_refpos.Popen')
    def test_unit_popencall( self, popen ):
        popen.return_value.stdout = 'tested'
        res = self._CM( self.bam, '', 25, 100000 )
        eq_( 'tested', res )
        popen.assert_called_with(['samtools','mpileup','-Q','25','-d','100000',self.bam],stdout=-1)
        self._CM( self.bam, 'den1:1-5', 25, 100000 )
        eq_( 'tested', res )
        popen.assert_called_with(['samtools','mpileup','-Q','25','-d','100000','-r','den1:1-5',self.bam],stdout=-1)

class TestStatsAtPos(Base):
    def setUp( self ):
        super(TestStatsAtPos,self).setUp()
        self.mp = {1046: join( fixtures.THIS, 'fixtures', 'mpileup_1046.txt' )}

    def _C( self, pos, bam, minmq=0, maxd=100000 ):
        from stats_at_refpos import stats_at_pos
        return stats_at_pos( pos, bam, minmq, maxd )

    @patch('stats_at_refpos.mpileup')
    def test_func_works( self, mpileup ):
        # Patch mpileup with expected stuffs
        mpileup.return_value = open(self.mp[1046])
        res = self._C( 1046, self.bam )
        eb = {
                'A':{'AvgReadQ':0.0,'AvgMapQ':36.02,'Depth':1349,'AvgReadQ':0.0,'PctTotal':25.49},
                'C':{'AvgReadQ':0.0,'AvgMapQ':14.25,'Depth':4,'AvgReadQ':0.0,'PctTotal':0.08},
                '*':{'AvgReadQ':0.0,'AvgMapQ':39.0,'Depth':1,'AvgReadQ':0.0,'PctTotal':0.02},
                'T':{'AvgReadQ':0.0,'AvgMapQ':17.5,'Depth':4,'AvgReadQ':0.0,'PctTotal':0.08},
                'G':{'AvgReadQ':0.0,'AvgMapQ':36.75,'Depth':3934,'AvgReadQ':0.0,'PctTotal':74.34}
            }
        e = {
            'Bases': eb,
            'AvgMapQ': 36.54,
            'AvgReadQ': 0.0,
            'TotalDepth': 5292
        }

        for base, valuesd in eb.iteritems():
            for k,v in valuesd.items():
                eq_( v, res['Bases'][base][k], "{} != {} for {}".format(v, res['Bases'][base][k], k) )

        for k in e:
            eq_( e[k], res[k] )
