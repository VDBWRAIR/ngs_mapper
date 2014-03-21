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
