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
        self.mp = {1046: join( fixtures.THIS, 'fixtures', 'mpileup_1046.txt' )}
        self.bam = self.__class__.bam

class TestMpileupPySam(Base):
    def _CM( self, bamfile, regionstr, minqual, maxd ):
        from stats_at_refpos import mpileup_pysam
        return mpileup_pysam( bamfile, regionstr, minqual, maxd )

    @patch('pysam.Samfile')
    def test_unit_call( self, samfile ):
        self._CM( self.bam, '', 25, 100000 )
        samfile.assert_called_with(self.bam)
        samfile.return_value.pileup.assert_called_with(region='')

class TestMpileupPopen(Base):
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

class StatsAtPos(Base):
    def _do( self, res ):
        from collections import OrderedDict
        eb = OrderedDict([
                ('G',{'AvgReadQ':0.0,'AvgMapQ':36.74,'Depth':3936,'AvgReadQ':0.0,'PctTotal':74.32}),
                ('A',{'AvgReadQ':0.0,'AvgMapQ':36.01,'Depth':1350,'AvgReadQ':0.0,'PctTotal':25.49}),
                ('C',{'AvgReadQ':0.0,'AvgMapQ':14.25,'Depth':4,'AvgReadQ':0.0,'PctTotal':0.08}),
                ('T',{'AvgReadQ':0.0,'AvgMapQ':17.5,'Depth':4,'AvgReadQ':0.0,'PctTotal':0.08}),
                ('*',{'AvgReadQ':0.0,'AvgMapQ':39.0,'Depth':1,'AvgReadQ':0.0,'PctTotal':0.02}),
                ('N',{'AvgReadQ':0.0,'AvgMapQ':2.0,'Depth':1,'AvgReadQ':0.0,'PctTotal':0.02})
            ])
        e = {
            'Bases': eb,
            'AvgMapQ': 36.52,
            'AvgReadQ': 0.0,
            'TotalDepth': 5296
        }

        # Ensure keys are same and in same order
        eq_( eb.keys(), res['Bases'].keys() )

        for base, valuesd in eb.iteritems():
            for k,v in valuesd.items():
                eq_( v, res['Bases'][base][k], "{} != {} for {}".format(v, res['Bases'][base][k], k) )

        for k in e:
            eq_( e[k], res[k] )

class TestStatsAtPosPysam(StatsAtPos):
    def _C( self, regionstr, bam, minmq=0, maxd=100000 ):
        from stats_at_refpos import stats_at_pos_pysam
        return stats_at_pos_pysam( regionstr, bam, minmq, maxd=100000 )

class TestStatsAtPosPopen(StatsAtPos):
    def _C( self, pos, bam, minmq=0, maxd=100000 ):
        from stats_at_refpos import stats_at_pos
        return stats_at_pos( pos, bam, minmq, maxd )

    @patch('stats_at_refpos.mpileup')
    def test_func_works( self, mpileup ):
        # Patch mpileup with expected stuffs
        mpileup.return_value = open(self.mp[1046])
        res = self._C( 1046, self.bam )
        self._do( res )
