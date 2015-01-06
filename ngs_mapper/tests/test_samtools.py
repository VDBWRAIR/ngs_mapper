import common
import fixtures

from nose.tools import eq_, raises, ok_
from nose.plugins.attrib import attr
from mock import MagicMock, patch, Mock, call

from os.path import *
import os

class Base(common.BaseBaseCaller):
    modulepath = 'ngs_mapper.samtools'

    def setUp( self ):
        super(Base,self).setUp()
        self.mp = {1046: join( fixtures.THIS, 'fixtures', 'mpileup_1046.txt' )}

class TestView(Base):
    functionname = 'view'
    '''
Usage:   samtools view [options] <in.bam>|<in.sam> [region1 [...]]

Options: -b       output BAM
         -h       print header for the SAM output
         -H       print header only (no alignments)
         -S       input is SAM
         -u       uncompressed BAM output (force -b)
         -1       fast compression (force -b)
         -x       output FLAG in HEX (samtools-C specific)
         -X       output FLAG in string (samtools-C specific)
         -c       print only the count of matching records
         -B       collapse the backward CIGAR operation
         -@ INT   number of BAM compression threads [0]
         -L FILE  output alignments overlapping the input BED FILE [null]
         -t FILE  list of reference names and lengths (force -S) [null]
         -T FILE  reference sequence file (force -S) [null]
         -o FILE  output file name [stdout]
         -R FILE  list of read groups to be outputted [null]
         -f INT   required flag, 0 for unset [0]
         -F INT   filtering flag, 0 for unset [0]
         -q INT   minimum mapping quality [0]
         -l STR   only output reads in library STR [null]
         -r STR   only output reads in read group STR [null]
         -s FLOAT fraction of templates to subsample; integer part as seed [-1]
         -?       longer help
    '''
    @patch('ngs_mapper.samtools.Popen')
    def test_returns_stdout( self, popen ):
        popen.return_value.stdout = 'tested'
        res = self._C( self.bam )
        eq_( 'tested', res )

    @patch('__builtin__.open')
    @patch('ngs_mapper.samtools.Popen')
    def test_inbam_outsam( self, popen, open ):
        open.return_value = 'null'
        res = self._C( self.bam )
        cmd = ['samtools','view',self.bam]
        popen.assert_called_with(cmd,stdout=-1)

    @patch('__builtin__.open')
    @patch('ngs_mapper.samtools.Popen')
    def test_pipeinput( self, popen, open ):
        open.return_value = 'null'
        infile = Mock()
        infile.fileno.return_value = -5
        res = self._C( infile )
        cmd = ['samtools','view','-']
        popen.assert_called_with(cmd,stdout=-1,stdin=-5)

    def test_bam_to_sam_and_back( self ):
        # Just pipe samtools sam output into samtools to bam output
        sam = self._C( self.bam, h=True, u=True )
        bam = self._C( sam, h=True, b=True )
        with open( 'new.bam', 'wb' ) as fh:
            fh.write( bam.read() )
        # Original file and new file should be same size
        eq_( os.stat(self.bam).st_size, os.stat('new.bam').st_size )

    def test_regionstring( self ):
        # Returns all reads for Ref2 since the test sam file 
        # has reads under every base(reads for ref 2 start at 10 and go to 20
        res = self._C( self.bam, 'Ref2:3-5' )
        i = 0
        for i, line in enumerate( res, 10 ):
            line = line.split()
            print line
            eq_( 'Read{}'.format(i), line[0] )
        eq_( 20, i )

class TestProp(Base):
    functionname = 'Prop'

    def test_does_type( self ):
        class A(object):
            a = self._C( 'i', int )
            b = self._C( 'f', float )
        a = A()
        a.a = '5'
        a.b = '5.1'
        eq_( 5, a.a )
        eq_( 5.1, a.b )

########### SamRow Tests ################
class SamRowBase(Base):
    functionname = 'SamRow'

    def setUp(self):
        super(SamRowBase,self).setUp()
        self.row = 'Read1	0	Ref1	1	60	8M	=	0	0	ACGTACGT	IIIIIIII	'

class TestSamRow(SamRowBase):
    def test_properties_set( self ):
        row = self.row + 'NM:i:0	AS:i:250'
        r = self._C( row )
        eq_( 'Read1', r.QNAME )
        eq_( 0, r.FLAG )
        eq_( 'Ref1', r.RNAME )
        eq_( 1, r.POS )
        eq_( 60, r.MAPQ )
        eq_( '8M', r.CIGAR )
        eq_( '=', r.RNEXT )
        eq_( 0, r.PNEXT )
        eq_( 0, r.TLEN )
        eq_( 'ACGTACGT', r.SEQ )
        eq_( 'IIIIIIII', r._qual )
        eq_( [40]*8, r.QUAL )
        eq_( [('NM',0),('AS',250)], r.TAGS )
        eq_( 'NM:i:0	AS:i:250', r._tags )

    def test_notags( self ):
        r = self._C( self.row[:-1] )
        eq_( '', r._tags )

class TestUnitSamRowStr(SamRowBase):
    def test_samestring_notags( self ):
        r = self._C( self.row[:-1] )
        eq_( self.row[:-1], str(r) )

    def test_samestring_tags( self ):
        r = self._C( self.row + 'NM:i:0' )
        eq_( self.row + 'NM:i:0', str(r) )

class TestUnitTagsToList(SamRowBase):
    def test_notags( self ):
        r = self._C( self.row )
        eq_( [], r.TAGS )
        eq_( '', r._tags )

    def test_single_and_convertstype( self ):
        types = (
            ('A',str),
            ('i',int),
            ('f',float),
            ('Z',str),
            ('H',hex),
            ('B',list)
        )
        for c,t in types:
            if c == 'B':
                row = self.row + 'aa:{}:'.format(c) + '1,1'
                r = self._C( row )
                eq_( ('aa',[1,1]), r.TAGS[0] )
            else:
                row = self.row + 'aa:{}:'.format(c) + str(t(1))
                r = self._C( row )
                if c == 'H':
                    t = str
                ok_( isinstance( r.TAGS[0][1], t ), "Expected {} to be a {} but got {}".format(c,t,type(r.TAGS[0][1])) )
            ok_( c in r._tags, "Did not set the correct field type" )

    def test_multi( self ):
        self.row += 'NM:i:0	AS:i:250'
        r = self._C( self.row )
        eq_( ('NM',0), r.TAGS[0] )
        eq_( ('AS',250), r.TAGS[1] )

class MpileupBase(Base):
    pass

class TestMpileup(MpileupBase):
    functionname = 'mpileup'

    @patch('__builtin__.open')
    @patch('ngs_mapper.samtools.Popen')
    def test_unit_popencall( self, popen, open ):
        open.return_value = 'null'
        popen.return_value.stdout = 'tested'
        res = self._C(self.bam, '', 20, 25, 100000 )
        eq_( 'tested', res )
        popen.assert_called_with(['samtools','mpileup','-s','-q','20','-Q','25','-d','100000',self.bam],stdout=-1, stderr='null')
        self._C( self.bam, 'den1:1-5', 20, 25, 100000 )
        eq_( 'tested', res )
        popen.assert_called_with(['samtools','mpileup','-s','-q','20','-Q','25','-d','100000','-r','den1:1-5',self.bam],stdout=-1, stderr='null')

    @attr('current')
    @patch('ngs_mapper.samtools.Popen')
    def test_called_with_individual_refs(self, mock_popen):
        mock_popen.side_effect = [Mock(stdout=self.mpileups[ref]) for ref in sorted(self.mpileups)]
        # Refname, reflen
        rets = []
        for ref in sorted(self.mpileups):
            rets.append(self._C(self.bam, ref, 0, 0, 100))
        for call, r in zip(mock_popen.call_args_list, rets):
            # This should be the reference name from the built cmd
            ref = call[0][0][10]
            eq_(self.mpileups[ref], r)

    @attr('current')
    @patch('ngs_mapper.samtools.Popen')
    def test_testbam_all_refs(self, mock_popen):
        _all = []
        for ref, piles in self.mpileups.iteritems():
            for pile in piles:
                _all.append(pile)
        mock_popen.return_value.stdout = _all

        r = self._C( self.bam, None, 0, 0, 100 )
        eq_(r, _all)

class TestNoGapMpileup(MpileupBase):
    functionname = 'nogap_mpileup'

    @attr('current')
    @patch('ngs_mapper.samtools.Popen')
    def test_mpileup_can_parse(self, mock_popen):
        mock_popen.return_value.stdout = self.mpileups['Ref5']
        from ngs_mapper.samtools import MPileupColumn
        for row in self._C(self.bam, 'Ref5:2-4', 0, 0, 100):
            MPileupColumn(row)

    @attr('current')
    @patch('ngs_mapper.samtools.Popen')
    def test_called_with_individual_refs(self, mock_popen):
        _all = [
            self._mock_pileup_str('R1',1,'A',1,'A','!','!'),
            self._mock_pileup_str('R1',3,'A',1,'A','!','!'),
            self._mock_pileup_str('R2',3,'A',1,'A','!','!'),
            self._mock_pileup_str('R2',5,'A',1,'A','!','!'),
        ]
        mock_popen.return_value.stdout = _all
        _ex = [
            self._mock_pileup_str('R1',1,'A',1,'A','!','!'),
            self._mock_pileup_str('R1',2,'',0,'','',''),
            self._mock_pileup_str('R1',3,'A',1,'A','!','!'),
            self._mock_pileup_str('R2',1,'',0,'','',''),
            self._mock_pileup_str('R2',2,'',0,'','',''),
            self._mock_pileup_str('R2',3,'A',1,'A','!','!'),
            self._mock_pileup_str('R2',4,'',0,'','',''),
            self._mock_pileup_str('R2',5,'A',1,'A','!','!'),
        ]
        r = self._C(self.bam, None, 0, 0, 100)
        for _e, _r in zip(_ex, list(r)):
            eq_(_e, _r)

    @attr('current')
    @patch('ngs_mapper.samtools.Popen')
    def test_fills_missing_positions(self, mock_popen):
        _all = [
            self._mock_pileup_str('R1', i, 'A', 1, 'A', '!', '!')
            for i in range(2, 11, 2)
        ]
        mock_popen.return_value.stdout = _all
        r = self._C(self.bam, None, 0, 0, 100)
        # Each ref should not have gaps
        expected = []
        for i in range(1, 11):
            if i % 2 == 0:
                expected.append(self._mock_pileup_str('R1', i, 'A', 1, 'A', '!', '!'))
            else:
                expected.append(self._mock_pileup_str('R1', i, '', 0, '', '', ''))

        for ex, re in zip(expected,list(r)):
            eq_(ex.split(), re.split())

class TestUnitCharToQual(Base):
    functionname = 'char_to_qual'

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
class MpileupBase(Base):
    functionname = 'MPileupColumn'

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

    def test_issue_174( self ):
        str = 'Den3/KDC0070A/Thailand/2010/Den3_1	88	N	36	*********A*********AAAA*********AAAt	CA?<=;=??;B?;;?;==<;;;D==<???=?;H;;:'
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

    def test_issue_174( self ):
        # Issue when there are no mapping qualities
        # mapq should be all set to 0
        str = 'Den3/KDC0070A/Thailand/2010/Den3_1	88	N	36	*********A*********AAAA*********AAAt	CA?<=;=??;B?;;?;==<;;;D==<???=?;H;;:'
        r = self._CA( str )
        eq_( [0]*8, r['A']['mapq'] )

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

class TestUnitMpileupStr(MpileupBase):
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

class TestUnitParseRegionString(Base):
    functionname = 'parse_regionstring'

    def test_invalid_regionstring(self):
        from ngs_mapper.samtools import InvalidRegionString
        bad_strings = ('', None, 'ref1:2-1')

        @raises(InvalidRegionString)
        def _test(string):
            self._C(string)

        for rstr in bad_strings:
            yield _test, rstr

    def test_start_gt_stop( self ):
        from ngs_mapper.samtools import InvalidRegionString
        try:
            self._C( 'ref1:2-1' )
            assert False, "Did not raise InvalidRegionString"
        except InvalidRegionString as e:
            assert True

    def test_incorrect_format( self ):
        from ngs_mapper.samtools import InvalidRegionString
        try:
            self._C( 'sometext' )
            self._C( 'sometext:1' )
            self._C( 'sometext:1-' )
            self._C( 'sometext:1- ' )
            self._C( 'sometext:a-b' )
            assert False, "Did not raise InvalidRegionString"
        except InvalidRegionString as e:
            assert True

    def test_correct_singlebase( self ):
        r = self._C( 'ref:1-1' )
        eq_( ('ref',1,1), r )

    def test_correct_multibase( self ):
        r = self._C( 'ref:1-2' )
        eq_( ('ref',1,2), r )
