import fixtures
import common
from bam import indexbam

from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import MagicMock, patch, Mock, call

from os.path import *
import os
from collections import OrderedDict

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
    @patch('stats_at_refpos.pysam_col')
    def test_unit_call( self, pysam_col, samfile ):
        samfile.return_value.pileup.return_value = [Mock(pos=1)]
        res = self._CM( self.bam, 'Ref1:1-1', 25, 100000 )
        next(res)
        samfile.assert_called_with(self.bam)
        samfile.return_value.pileup.assert_called_with(region='Ref1:1-1')

    def test_func_call( self ):
        res = self._CM( self.bam, 'Den1/U88535_1/WestPac/1997/Den1_1:7935-7935', 0, 10000 )
        e = 'Den1/U88535_1/WestPac/1997/Den1_1\t7935\tN\t19\tGGGGGGGGGGGGGGGGGGG\tFGFBHCFHHHHH1HHFBHF\t<<<<<<<<<<<<<<<<<<<'
        eq_( e, next(res) )

class TestPysamCol(Base):
    def _C( self, pysamcol, refname ):
        from stats_at_refpos import pysam_col
        return pysam_col( pysamcol, refname )

    def mock_pileup_read( self, qpos, mapq, quals, seq, pos ):
        return Mock(
            qpos = qpos,
            alignment = Mock(
                mapq = mapq,
                qual = quals,
                seq=seq,
                pos = pos
            )
        )

    def test_func_pysamcol( self ):
        # Mock spot 10 on mock ref
        ref = 'Ref1'
        # 2 A, 8 G
        col = [
            self.mock_pileup_read( 2, '<', 'I'*10, 'A'*8+'G'*2, 8 ),
            self.mock_pileup_read( 2, '<', 'I'*10, 'A'*8+'G'*2, 8 ),
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, '<', 'I'*10, 'A'*8+'G'*2, 1 ), 
        ]
        res = self._C( col, ref )
        eq_( 
            '\t'.join( [ref,'10','N','10','AAGGGGGGGG','I'*10,'<'*10] ),
            res
        )

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
        eb = OrderedDict([
                ('G',{'AvgReadQ':0.0,'AvgMapQ':36.74,'Depth':3936,'PctTotal':74.32}),
                ('A',{'AvgReadQ':0.0,'AvgMapQ':36.01,'Depth':1350,'PctTotal':25.49}),
                ('C',{'AvgReadQ':0.0,'AvgMapQ':14.25,'Depth':4,'PctTotal':0.08}),
                ('T',{'AvgReadQ':0.0,'AvgMapQ':17.5,'Depth':4,'PctTotal':0.08}),
                ('*',{'AvgReadQ':0.0,'AvgMapQ':39.0,'Depth':1,'PctTotal':0.02}),
                ('N',{'AvgReadQ':0.0,'AvgMapQ':2.0,'Depth':1,'PctTotal':0.02})
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

class TestCompileStats(Base):
    def _C( self, stats ):
        from stats_at_refpos import compile_stats
        return compile_stats( stats )

    def test_func_works( self ):
        stats = {
            'depth': 1000,
            'mqualsum': 30*900+40*100,
            'rqualsum': 0,
            'G': [30]*900,
            'A': [40]*100
        }
        res = self._C( stats )

        eq_( stats['depth'], res['TotalDepth'] )
        eq_( 31.0, res['AvgMapQ'] )
        eq_( 0.0, res['AvgReadQ'] )

        g = res['Bases']['G']
        eq_( 900, g['Depth'] )
        eq_( 30.0, g['AvgMapQ'] )
        eq_( 0.0, g['AvgReadQ'] )
        eq_( 90.0, g['PctTotal'] )

        a = res['Bases']['A']
        eq_( 100, a['Depth'] )
        eq_( 40.0, a['AvgMapQ'] )
        eq_( 0.0, a['AvgReadQ'] )
        eq_( 10.0, a['PctTotal'] )
