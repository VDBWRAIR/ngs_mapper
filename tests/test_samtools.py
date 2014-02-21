import common
import fixtures

from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import MagicMock, patch, Mock, call

from os.path import *
import os

class Base(common.BaseBamRef):
    def setUp( self ):
        super(Base,self).setUp()
        self.mp = {1046: join( fixtures.THIS, 'fixtures', 'mpileup_1046.txt' )}
        self.bam = self.__class__.bam

class TestMpileup(Base):
    def _CM( self, bamfile, regionstr, minmq, minbq, maxd ):
        from samtools import mpileup
        return mpileup( bamfile, regionstr, minmq, minbq, maxd )

    @patch('__builtin__.open')
    @patch('samtools.Popen')
    def test_unit_popencall( self, popen, open ):
        open.return_value = 'null'
        popen.return_value.stdout = 'tested'
        res = self._CM( self.bam, '', 20, 25, 100000 )
        eq_( 'tested', res )
        popen.assert_called_with(['samtools','mpileup','-s','-q','20','-Q','25','-d','100000',self.bam],stdout=-1, stderr='null')
        self._CM( self.bam, 'den1:1-5', 20, 25, 100000 )
        eq_( 'tested', res )
        popen.assert_called_with(['samtools','mpileup','-s','-q','20','-Q','25','-d','100000','-r','den1:1-5',self.bam],stdout=-1, stderr='null')


class TestUnitCharToQual(object):
    def _C( self, qual_char ):
        from samtools import char_to_qual
        return char_to_qual( qual_char )

    def test_correct_phred( self ):
        # Ascii ! -> 33, ] -> 93
        # Check all qualities from 0 -> 60
        for i in range( 0, 60+1 ):
            # Phred is + 33
            c = chr( i + 33 )
            e = ord( c )
            yield self.check_correct, i, c

    def check_correct( self, e, c ):
        eq_( e, self._C( c ) )

########### MPileupColumn Tests ################
class MpileupBase(common.BaseClass):
    def _C( self, pileupstring ):
        from samtools import MPileupColumn
        return MPileupColumn( pileupstring )

class TestUnitBases(MpileupBase):
    def test_lowercaseuppercase( self ):
        str = 'Ref1	1	N	10	AaCcGgTtNn	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 'AACCGGTTNN', r.bases )

    def test_dotcomma( self ):
        str = 'Ref1	1	A	12	.,CcGgTtNn.,	IIIIIIIIIIII	]]]]]]]]]]]]'
        r = self._C( str )
        eq_( 'AACCGGTTNNAA', r.bases )

    def test_insertdelete( self ):
        str = 'Ref1	1	A	10	G+2AA-2AA*.G,.....	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 'G*AGAAAAAA', r.bases )

    def test_endreadbeginread( self ):
        str = 'Ref1	1	N	10	A^]A$AAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 'AAAAAAAAAA', r.bases )

class TestUnitBQuals(MpileupBase):
    pass

class TestUnitMQuals(MpileupBase):
    def test_mquals_eq_bquals_len( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( [60]*10, r.mquals )
    
    def test_mquals_ne_bquals_len_same( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]]'
        r = self._C( str )
        eq_( [60]*10, r.mquals )

    def test_mquals_ne_bquals_len_diff( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]I]]]]]'
        r = self._C( str )
        eq_( [], r.mquals )

    def test_mquals_missing( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII'
        r = self._C( str )
        eq_( [], r.mquals )

class TestUnitIter(MpileupBase):
    def test_iterates( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        c = 0
        for i in r:
            eq_( i, ('A',40,60) )
            c += 1
        eq_( 10, c )

    def test_iterates_missingmquals( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII'
        r = self._C( str )
        c = 0
        for i in r:
            eq_( i, ('A',40,0) )
            c += 1
        eq_( 10, c )

class TestUnitBaseStats(MpileupBase):
    def _CA( self, mpstr ):
        return self._C( mpstr ).base_stats()

    def test_qualsums_set( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._CA( str )
        eq_( 40*10.0, r['bqualsum'] )
        eq_( 60*10.0, r['mqualsum'] )
        eq_( float, type(r['bqualsum']) )
        eq_( float, type(r['mqualsum']) )

    def test_depth_set( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._CA( str )
        eq_( 10, r['depth'] )

    def test_bases_set( self ):
        str = 'Ref1	1	A	14	Aa.,CcGgTtNn*$C^]	EDCBAIHGFEDCBA	ABCDEFGHIABCDE'
        r = self._CA( str )
        eq_( r['A']['baseq'], [36,35,34,33] )
        eq_( r['A']['mapq'], [32,33,34,35] )

class TestUnitAvgQuals(MpileupBase):
    def test_avgbqual_set( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	ABCDEABCDE	]]]]]]]]]]'
        r = self._C( str )
        eq_( 34.0, r.bqual_avg() )

    def test_example_1( self ):
        str =   'Den1/U88535_1/WestPac/1997/Den1_1	6109	N	13	GgnGgggggtGgg	CB#GHHHHG2GHH'
        r = self._C( str )
        eq_( 33.38, r.bqual_avg() )

    def test_avgmqual_set( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	ABCDEABCDE'
        r = self._C( str )
        eq_( 34.0, r.mqual_avg() )

class TestUnitStr(MpileupBase):
    def test_returns_correct_string( self ):
        s = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	ABCDEABCDE'
        r = self._C( s )
        eq_( str(r), s )

class TestUnitPileupColumnInit(MpileupBase):
    def test_sets_ref_pos_base_depth( self ):
        # 10 A, 40 BQ, 60 MQ
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 'Ref1', r.ref )
        eq_( 1, r.pos )
        eq_( 'N', r.refbase )
        eq_( 10, r.depth )
        eq_( 'A'*10, r._bases )
        eq_( 'I'*10, r._bquals )
        eq_( ']'*10, r._mquals )

    def test_old_samtools_mapqual_same( self ):
        # Older samtools did not trim out mapping qualities when -Q and -q were specified
        # If they were all the same then can just trim them down
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]]]]]]]]]]]'
        r = self._C( str )
        eq_( ']'*20, r._mquals )
        eq_( [60]*10, r.mquals )

    def test_old_samtools_mapqual_diff( self ):
        # otherwise they cannot be used so it should be set to None
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII	]]]]]A]]]]]]]]]]]]]]'
        r = self._C( str )
        eq_( ']]]]]A]]]]]]]]]]]]]]', r._mquals )
        eq_( [], r.mquals )

    def test_no_mapping_qualities( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA	IIIIIIIIII'
        r = self._C( str )
        eq_( [], r.mquals )
        eq_( 'A'*10, r.bases )
        eq_( [40]*10, r.bquals )

    def test_insert_start( self ):
        str = 'Ref1	1	N	10	+2NNAAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'A'*10, r.bases )

    def test_insert_end( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA+2NN	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'A'*10, r.bases )

    def test_insert_middle( self ):
        str = 'Ref1	1	N	10	AAAAA+2NNAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'A'*10, r.bases )

    def test_delete_start( self ):
        str = 'Ref1	1	N	10	-2NNAAAAAAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'A'*10, r.bases )

    def test_delete_end( self ):
        str = 'Ref1	1	N	10	AAAAAAAAAA-2NN	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'A'*10, r.bases )

    def test_delete_middle( self ):
        str = 'Ref1	1	N	10	AAAAA-2NNAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'A'*10, r.bases )

    def test_other_characters( self ):
        str = 'Ref1	1	N	10	A$.,aAA^]AA*A	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'I'*10, r._bquals )
        eq_( ']'*10, r._mquals )
        eq_( 'ANNAAAAA*A', r.bases )

    def test_indel_gt_9( self ):
        str = 'Ref1	1	N	10	AAAAA-10NNNNNNNNNNAAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'A'*10, r.bases )

    def test_gap( self ):
        # Guess gaps don't have quality scores
        str = 'Ref1	1	N	10	AAAAA-AAAAA	IIIIIIIIII	]]]]]]]]]]'
        r = self._C( str )
        eq_( 10, r.depth )
        eq_( 'I'*10, r._bquals )
        eq_( ']'*10, r._mquals )
        eq_( 'A'*10, r.bases )
