import common
from fixtures import FIXDIR

from mock import Mock, MagicMock, patch, call
from nose.tools import eq_, raises, timed
from nose.plugins.attrib import attr

from os.path import *
from StringIO import StringIO
import sys
import json

# For patching because I don't know a better way
PKG_BASE = __package__.rsplit('.', 1)[0]

t = Mock()
t.LCDT = 10
t.LCAQT = 10
t.GDT = 5
t.GAQT = 5

class Base(common.BaseTester):
    def setUp( self ):
        self.samplename = '00103-01'
        self.bamfile = join(FIXDIR,'jsonfix','00103-01.bam')
        self.jsonfile = join(FIXDIR,'jsonfix','00103-01.json')
        self.refstats = {
            'Den4/AY618992_1/Thailand/2001/Den4_1':
                ['Den4/AY618992_1/Thailand/2001/Den4_1','10649','147751','220'],
            '*':
                ['*','0','0','23842']
        }

class TestIntegrate(Base):
    def setUp(self):
        s = super(TestIntegrate,self)
        if hasattr(s,'setUp'):
            s.setUp()

    def test_get_refstats(self):
        from ngs_mapper.bam import get_refstats
        rs = get_refstats( self.bamfile )
        eq_( self.refstats, rs, "{} != {}. You may want to try to run nose with -s".format(self.refstats,rs) )

    def _cmpcontents(self, file1, file2 ):
        from difflib import unified_diff
        f1 = file1.read()
        f2 = file2.getvalue()+'\n'
        assert f2, 'Result file was empty'
        if f1.strip() == f2.strip():
            assert True
 
        dif = unified_diff( f1.splitlines(1), f2.splitlines(1), fromfile='expected', tofile='result' )
        s = ''.join(dif) + '\n' 
        if s.rstrip():
            sys.stderr.write( s )

        assert len(f1) > 0
        eq_( len(f1), len(f2) )
        eq_( f1, f2 )

    def test_parse_mapstats(self):
        from BamCoverage.bam import alignment_info
        from BamCoverage.bam_to_json import parse_mapstats
        import json
        regions = json.load(open(self.jsonfile))
        regions = regions['references']['Den4/AY618992_1/Thailand/2001/Den4_1']['regions']
        exp = [(r[0],r[1],str(r[2])) for r in regions]
        ai = alignment_info( self.bamfile )
        res = parse_mapstats( next(ai)[1] )
        for e,r in zip(exp,res):
            if e != r:
                print exp
                print res
            eq_( e, r )
        # Gap at end?
        if exp != res and exp[:-1] != res:
            eq_( exp, res )

    @patch('ngs_mapper.bam.get_refstats')
    @attr('slow')
    def test_output_json(self,refstats):
        refstats.return_value = self.refstats
        from BamCoverage.bam import alignment_info
        from BamCoverage.bam_to_json import output_json
        ai = alignment_info( self.bamfile )
        with patch('sys.stdout',new_callable=StringIO) as stdout:
            output_json( ai, self.samplename )
            self._cmpcontents( open(self.jsonfile), stdout )

# Integratish testing
@patch('BamCoverage.bam_to_json.Thresholds',t)
class TestParseMapstats(common.BaseTester):
    def _CPM( self, stats ):
        from BamCoverage.bam_to_json import parse_mapstats
        return parse_mapstats( stats )
    
    def test_beginning_gap(self):
        stats = {'depths':[0,0,0,100,100],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'Gap'),(4,5,'Normal')], self._CPM( stats ) )

    def test_beginning_lc(self):
        stats = {'depths':[6,6,6,100,100],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'LowCoverage'),(4,5,'Normal')], self._CPM( stats ) )

    def test_beginning_lq(self):
        stats = {'depths':[100,100,100,100,100],'avgquals':[6,6,6,40,40]}
        eq_( [(1,4,'LowQuality'),(4,5,'Normal')], self._CPM( stats ) )

    def test_beginning_n(self):
        stats = {'depths':[100,100,100,0,0],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'Normal'),(4,5,'Gap')], self._CPM( stats ) )


    def test_middle_gap(self):
        stats = {'depths':[100,0,0,100,100],'avgquals':[40,40,40,40,40]}
        eq_( [(1,2,'Normal'),(2,4,'Gap'),(4,5,'Normal')], self._CPM( stats ) )

    def test_middle_lc(self):
        stats = {'depths':[100,6,6,100,100],'avgquals':[40,40,40,40,40]}
        eq_( [(1,2,'Normal'),(2,4,'LowCoverage'),(4,5,'Normal')], self._CPM( stats ) )

    def test_middle_lq(self):
        stats = {'depths':[100,100,100,100,100],'avgquals':[40,6,6,40,40]}
        eq_( [(1,2,'Normal'),(2,4,'LowQuality'),(4,5,'Normal')], self._CPM( stats ) )

    def test_middle_n(self):
        stats = {'depths':[0,100,100,0,0],'avgquals':[40,40,40,40,40]}
        eq_( [(1,2,'Gap'),(2,4,'Normal'),(4,5,'Gap')], self._CPM( stats ) )


    def test_end_gap(self):
        stats = {'depths':[100,100,100,0,0],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'Normal'),(4,5,'Gap')], self._CPM( stats ) )

    def test_end_lc(self):
        stats = {'depths':[100,100,100,6,6],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'Normal'),(4,5,'LowCoverage')], self._CPM( stats ) )

    def test_end_lq(self):
        stats = {'depths':[100,100,100,100,100],'avgquals':[40,40,40,6,6]}
        eq_( [(1,4,'Normal'),(4,5,'LowQuality')], self._CPM( stats ) )

    def test_end_n(self):
        stats = {'depths':[0,0,0,100,100],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'Gap'),(4,5,'Normal')], self._CPM( stats ) )


    def test_g_to_lc(self):
        stats = {'depths':[0,0,0,6,6],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'Gap'),(4,5,'LowCoverage')], self._CPM( stats ) )

    def test_g_to_lq(self):
        stats = {'depths':[0,0,0,100,100],'avgquals':[40,40,6,6,6]}
        eq_( [(1,4,'Gap'),(4,5,'LowQuality')], self._CPM( stats ) )

    def test_g_to_n(self):
        pass

    def test_lc_to_g(self):
        stats = {'depths':[6,6,6,0,0],'avgquals':[40,40,40,40,40]}
        eq_( [(1,4,'LowCoverage'),(4,5,'Gap')], self._CPM( stats ) )

    def test_lc_to_lq(self):
        stats = {'depths':[6,6,6,100,100],'avgquals':[40,40,40,6,6]}
        eq_( [(1,4,'LowCoverage'),(4,5,'LowQuality')], self._CPM( stats ) )

    def test_lc_to_n(self):
        pass

    def test_lq_to_g(self):
        stats = {'depths':[100,100,0,0,0],'avgquals':[6,6,40,40,40]}
        eq_( [(1,3,'LowQuality'),(3,5,'Gap')], self._CPM( stats ) )
        
    def test_lq_to_lc(self):
        stats = {'depths':[100,100,100,6,6],'avgquals':[6,6,6,40,40]}
        eq_( [(1,4,'LowQuality'),(4,5,'LowCoverage')], self._CPM( stats ) )
    
    def test_lq_to_n(self):
        pass

    def test_n_to_g(self):
        pass
    
    def test_n_to_lc(self):
        pass
    
    def test_n_to_lq(self):
        pass

    def test_lc_trumps_lq(self):
        stats = {'depths':[100,100,6,6,100],'avgquals':[40,40,6,6,40]}
        eq_( [(1,3,'Normal'),(3,5,'LowCoverage'),(5,5,'Normal')], self._CPM( stats ) )


class TestGetRefstats(common.BaseTester):
    def _deq( self, d1, d2 ):
        # ensure same reference keys
        eq_( sorted(d1.keys()), sorted(d2.keys()) )
        # Ensure refstats are equal
        for d1r, d2r in zip(d1.keys(),d2.keys()):
            d1list = sorted( d1[d1r] )
            d2list = sorted( d2[d2r] )
            eq_( d1list, d2list )

    def _CGR( self, bam ):
        from ngs_mapper.bam import get_refstats
        return get_refstats( bam )

    @patch('ngs_mapper.bam.Popen')
    def test_parses( self, popen ):
        popen.return_value.communicate.return_value = (
            'ref1\t5\t100\t10\n' \
            'ref2\t10\t100\t10\n' \
            '*\t0\t0\t0\n',
            ''
        )
        res = self._CGR( 'test' )
        exp = {
            'ref1': ['ref1','5','100','10'],
            'ref2': ['ref2','10','100','10'],
            '*': ['*','0','0','0']
        }
        self._deq( exp, res )

class TestParsePileup(common.BaseTester):
    def _deq( self, d1, d2 ):
        # ensure same reference keys
        eq_( sorted(d1.keys()), sorted(d2.keys()) )
        # Ensure each reference stats are same
        for d1r, d2r in zip(d1.keys(),d2.keys()):
            d1items = sorted(d1[d1r].items())
            d2items = sorted(d2[d2r].items())
            eq_( d1items, d2items )

    def _CPP( self, pileup ):
        from ngs_mapper.bqd import parse_pileup
        return parse_pileup( pileup )

    def test_missingpositions(self):
        exp = {
            'ref1': dict(
                mind=0,maxd=3,minq=0.0,maxq=34.0,length=9,depths=[0,0,0,0,0,1,0,0,3],avgquals=[0.0,0.0,0.0,0.0,0.0,32.0,0.0,0.0,33.0]
            )
        }
        lines = [
            'ref1\t6\tN\t1\tA\tA\n',
            'ref1\t9\tN\t3\tAAA\tABC\n',
        ]
        res = self._CPP( lines )
        self._deq( exp, res )

    def test_singleref(self):
        exp = {
            'ref1': dict(
                mind=1,maxd=5,minq=32.0,maxq=36.0,length=4,depths=[1,5,3,3],avgquals=[32.0,34.0,33.0,33.0]
            )
        }
        lines = [
            'ref1\t1\tN\t1\tA\tA\n',
            'ref1\t2\tN\t5\tAAAAA\tABCDE\n',
            'ref1\t3\tN\t3\tAAA\tABC\n',
            'ref1\t4\tN\t3\tAAA\tABC\n',
        ]
        res = self._CPP( lines )
        self._deq( exp, res )

    def test_multiref(self):
        exp = {
            'ref1': dict(
                mind=1,maxd=5,minq=32.0,maxq=36.0,length=4,depths=[1,5,3,3],avgquals=[32.0,34.0,33.0,33.0]
            ),
            'ref2': dict(
                mind=1,maxd=1,minq=32.0,maxq=32.0,length=1,depths=[1],avgquals=[32.0]
            )
        }
        lines = [
            'ref1\t1\tN\t1\tA\tA\n',
            'ref1\t2\tN\t5\tAAAAA\tABCDE\n',
            'ref1\t3\tN\t3\tAAA\tABC\n',
            'ref1\t4\tN\t3\tAAA\tABC\n',
            'ref2\t1\tN\t1\tA\tA\n',
        ]
        res = self._CPP( lines )
        self._deq( exp, res )

class TestMpileup(common.BaseTester):
    cmd = ['samtools','mpileup','-d','100000']

    def _RMP( self, bam, rg ):
        from ngs_mapper.bqd import mpileup
        return mpileup( bam, rg )

    @patch('ngs_mapper.bqd.Popen')
    def test_retvalue(self,popen):
        popen.return_value.stdout = 'great'
        res = self._RMP( 'test', None )
        eq_( res, 'great' )

    @patch('ngs_mapper.bqd.Popen')
    def test_runs_mpileup_nors(self,popen):
        res = self._RMP( 'test', None )
        ecmd = self.cmd + ['test']
        eq_( [call(ecmd, stdout=-1)], popen.call_args_list )

    @patch('ngs_mapper.bqd.Popen')
    def test_runs_mpileup_rs(self,popen):
        res = self._RMP( 'test', 'test2' )
        ecmd = self.cmd + ['-r','test2','test']
        eq_( [call(ecmd, stdout=-1)], popen.call_args_list )
