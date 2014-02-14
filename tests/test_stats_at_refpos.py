import common
import fixtures

from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import MagicMock, patch, Mock, call

from os.path import *
import os
from collections import OrderedDict

class Base(common.BaseBamRef):
    def setUp( self ):
        super(Base,self).setUp()
        self.mp = {1046: join( fixtures.THIS, 'fixtures', 'mpileup_1046.txt' )}
        self.bam = self.__class__.bam

class TestMpileupPySam(Base):
    def _CM( self, bamfile, regionstr, minmq, minbq, maxd ):
        from stats_at_refpos import mpileup_pysam
        return mpileup_pysam( bamfile, regionstr, minmq, minbq, maxd )

    @patch('pysam.Samfile')
    @patch('stats_at_refpos.pysam_col')
    def test_unit_call( self, pysam_col, samfile ):
        samfile.return_value.pileup.return_value = [Mock(pos=1)]
        res = self._CM( self.bam, 'Ref1:1-1', 25, 0, 100000 )
        next(res)
        samfile.assert_called_with(self.bam)
        samfile.return_value.pileup.assert_called_with(region='Ref1:1-1')

    def test_func_call( self ):
        res = self._CM( self.bam, 'Den1/U88535_1/WestPac/1997/Den1_1:7935-7935', 0, 0, 10000 )
        e = 'Den1/U88535_1/WestPac/1997/Den1_1\t7935\tN\t19\tGGGGGGGGGGGGGGGGGGG\tFGFBHCFHHHHH1HHFBHF\t<<<<<<<<<<<<<<<<<<<'
        eq_( e, next(res) )

class TestPysamCol(Base):
    def _C( self, pysamcol, refname, minmq, minbq ):
        from stats_at_refpos import pysam_col
        return pysam_col( pysamcol, refname, minmq, minbq )

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
    
    def test_missing_mq( self ):
        # Sometimes the mapping quality is missing for some reason
        ref = 'Ref1'
        col = [
            self.mock_pileup_read( 2, None, 'I'*10, 'A'*10, 8 ),
            self.mock_pileup_read( 2, '', 'I'*10, 'A'*10, 8 ),
            self.mock_pileup_read( 2, 32, 'I'*10, 'A'*10, 8 ),
            self.mock_pileup_read( 2, ']', 'I'*10, 'A'*10, 8 ),
            self.mock_pileup_read( 2, ']', 'I'*10, 'A'*10, 8 ),
        ]
        r = self._C( col, ref, 0, 0 )
        e = '\t'.join( [ref,'10','N','5','A'*5,'I'*5,'!!!]]'] )
        eq_( e, r )

    def test_func_allreads( self ):
        # Mock spot 10 on mock ref
        ref = 'Ref1'
        # 2 A, 8 G
        # ord('<') == 60
        # ord('I') == 73
        # ord('!') == 33
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
        res = self._C( col, ref, 0, 0 )
        eq_( 
            '\t'.join( [ref,'10','N','10','AAGGGGGGGG','I'*10,'<'*10] ),
            res
        )

    def test_func_minmq( self ):
        ref = 'Ref1'
        minmq = 25
        keepmq = chr(minmq) # Should be a :
        filtermq = chr(minmq-1)
        col = [
            self.mock_pileup_read( 2, keepmq, 'I'*10, 'A'*8+'G'*2, 8 ),
            self.mock_pileup_read( 2, keepmq, 'I'*10, 'A'*8+'G'*2, 8 ),
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, filtermq, 'I'*10, 'A'*8+'G'*2, 1 ), 
        ]
        res = self._C( col, ref, keepmq, 0 )
        eq_( 
            '\t'.join( [ref,'10','N','2','AA','I'*2,keepmq*2] ),
            res
        )

    def test_func_minbq( self ):
        ref = 'Ref1'
        minbq = 25
        keep = chr(minbq+33)
        filter = chr(minbq+33-1)
        mq = '<'
        col = [
            self.mock_pileup_read( 2, mq, keep*10, 'A'*8+'G'*2, 8 ),
            self.mock_pileup_read( 2, mq, keep*10, 'A'*8+'G'*2, 8 ),
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
            self.mock_pileup_read( 9, mq, filter*10, 'A'*8+'G'*2, 1 ), 
        ]
        res = self._C( col, ref, 0, minbq )
        eq_( 
            '\t'.join( [ref,'10','N','2','AA',keep*2,mq*2] ),
            res
        )

class TestMpileupPopen(Base):
    def _CM( self, bamfile, regionstr, minqual, maxd ):
        from stats_at_refpos import mpileup_popen
        return mpileup_popen( bamfile, regionstr, minqual, maxd )

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
    def _doit( self, res, eb, e ):
        # Ensure keys are same and in same order
        eq_( eb.keys(), res['Bases'].keys() )

        for base, valuesd in eb.iteritems():
            for k,v in valuesd.items():
                eq_( v, res['Bases'][base][k], "{} != {} for {}".format(v, res['Bases'][base][k], k) )

        for k in e:
            eq_( e[k], res[k] )

class TestStatsAtPos(StatsAtPos):
    def _C( self, bam, regionstr, minmq=0, minbq=0, maxd=100000 ):
        from stats_at_refpos import stats_at_pos
        return stats_at_pos( bam, regionstr, minmq, minbq, maxd )

    def test_func_works( self ):
        ref = 'Den1/U88535_1/WestPac/1997/Den1_1'
        pos = '6109'
        regionstr = '{}:{}-{}'.format(ref,pos,pos)
        res = self._C( self.bam, regionstr, 0, 0, 100000 )
        #Den1/U88535_1/WestPac/1997/Den1_1  6109    N   13  GgnGgggggtGgg   CB#GHHHHG2GHH
        eb = OrderedDict([
                ('G',{'AvgBaseQ':37.73,'AvgMapQ':60.0,'Depth':11,'PctTotal':84.62}),
                ('T',{'AvgBaseQ':17.0,'AvgMapQ':60.0,'Depth':1,'PctTotal':7.69}),
                ('N',{'AvgBaseQ':2.0,'AvgMapQ':60.0,'Depth':1,'PctTotal':7.69})
            ])
        e = {
            'Bases': eb,
            'AvgMapQ': 60.0,
            'AvgBaseQ': 33.38,
            'TotalDepth': 13
        }
        self._doit( res, eb, e )

class TestCompileStats(Base):
    def _C( self, stats ):
        from stats_at_refpos import compile_stats
        return compile_stats( stats )

    def test_func_works( self ):
        stats = {
            'depth': 1000,
            'mqualsum': 50*900+60*100,
            'bqualsum': 30*900+40*100,
            'G': {'mapq': [50]*900, 'baseq': [30]*900},
            'A': {'mapq': [60]*100, 'baseq': [40]*100}
        }
        res = self._C( stats )

        eq_( stats['depth'], res['TotalDepth'] )
        eq_( 51.0, res['AvgMapQ'] )
        eq_( 31.0, res['AvgBaseQ'] )

        g = res['Bases']['G']
        eq_( 900, g['Depth'] )
        eq_( 50.0, g['AvgMapQ'] )
        eq_( 30.0, g['AvgBaseQ'] )
        eq_( 90.0, g['PctTotal'] )

        a = res['Bases']['A']
        eq_( 100, a['Depth'] )
        eq_( 60.0, a['AvgMapQ'] )
        eq_( 40.0, a['AvgBaseQ'] )
        eq_( 10.0, a['PctTotal'] )

class TestMain(StatsAtPos):
    def _CM( self, args ):
        from stats_at_refpos import main
        return main( args )

    @patch('stats_at_refpos.stats_at_pos')
    def test_unit_runs( self, stats_at_pos ):
        args = Mock(
            regionstr='ref1:1-1',
            bamfile='somefile.bam',
            minmq=0,
            minbq=0,
            maxd=100000
        )
        self._CM( args )

    def test_func_runs( self ):
        args = Mock(
            bamfile=self.bam,
            regionstr='Den1/U88535_1/WestPac/1997/Den1_1:6109-6109',
            minmq=0,
            minbq=0,
            maxd=100000
        )
        eb = OrderedDict([
                ('G',{'AvgBaseQ':37.73,'AvgMapQ':60.0,'Depth':11,'PctTotal':84.62}),
                ('T',{'AvgBaseQ':17.0,'AvgMapQ':60.0,'Depth':1,'PctTotal':7.69}),
                ('N',{'AvgBaseQ':2.0,'AvgMapQ':60.0,'Depth':1,'PctTotal':7.69})
            ])
        e = {
            'Bases': eb,
            'AvgMapQ': 60.0,
            'AvgBaseQ': 33.38,
            'TotalDepth': 13
        }
        res = self._CM( args )
        self._doit( res, eb, e )

    def test_func_filters_minmq( self ):
        args = Mock(
            bamfile=self.bam,
            regionstr='Den1/U88535_1/WestPac/1997/Den1_1:6109-6109',
            minmq=61,
            minbq=0,
            maxd=100000
        )
        res = self._CM( args )
        #Den1/U88535_1/WestPac/1997/Den1_1  6109    N   13  GgnGgggggtGgg   CB#GHHHHG2GHH
        eb = OrderedDict([
                #('G',{'AvgBaseQ':0.0,'AvgMapQ':37.73,'Depth':11,'PctTotal':84.62}),
                #('T',{'AvgBaseQ':0.0,'AvgMapQ':17.0,'Depth':1,'PctTotal':7.69}),
                #('N',{'AvgBaseQ':0.0,'AvgMapQ':2.0,'Depth':1,'PctTotal':7.69})
            ])
        e = {
            'Bases': eb,
            'AvgMapQ': 0.0,
            'AvgBaseQ': 0.0,
            'TotalDepth': 0
        }
        self._doit( res, eb, e )

    def test_func_filters_minbq( self ):
        args = Mock(
            bamfile=self.bam,
            regionstr='Den1/U88535_1/WestPac/1997/Den1_1:6109-6109',
            minmq=0,
            minbq=30,
            maxd=100000
        )
        res = self._CM( args )
        #Den1/U88535_1/WestPac/1997/Den1_1  6109    N   13  GgnGgggggtGgg   CB#GHHHHG2GHH
        eb = OrderedDict([
                ('G',{'AvgBaseQ':37.73,'AvgMapQ':60.0,'Depth':11,'PctTotal':100.0}),
                #('T',{'AvgBaseQ':0.0,'AvgMapQ':17.0,'Depth':1,'PctTotal':7.69}),
                #('N',{'AvgBaseQ':0.0,'AvgMapQ':2.0,'Depth':1,'PctTotal':7.69})
            ])
        e = {
            'Bases': eb,
            'AvgMapQ': 60.0,
            'AvgBaseQ': 37.73,
            'TotalDepth': 11
        }
        self._doit( res, eb, e )
