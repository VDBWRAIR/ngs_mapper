import common
import fixtures

from nose.tools import eq_, raises, timed
from nose.plugins.attrib import attr
from mock import patch, Mock, MagicMock

from StringIO import StringIO
from os.path import join, dirname, basename, exists

# How long you expect each base position on the reference to take to process
EXPECTED_TIME_PER_BASE = 0.0015
# How long you expect for the entire fixtures/base_caller/test.bam to complete
EXPECTED_TOTAL_TIME = 0.6

class Base( common.BaseBaseCaller ):
    def mock_stats( self ):
        s = { 'baseq': [40]*10, 'mapq': [60]*10 }
        self.stats = self.make_stats({
            b: s.copy() for b in 'ATGC'
        })

    def make_stats( self, base_stats ):
        stats = {}
        for k,v in base_stats.items():
            stats[k] = {}
            stats[k]['baseq'] = v['baseq']
            if 'mapq' not in stats[k]:
                stats[k]['mapq'] = v['baseq']
            else:
                stats[k]['mapq'] = v['mapq']

        self.update_stats( stats )

        return stats

    def update_stats( self, stats ):
        ''' Updates depth, mqualsum and bqualsum from a given stats '''
        nonbasekeys = ('depth','mqualsum','bqualsum')
        # Reset to 0
        for k in nonbasekeys:
            stats[k] = 0
        
        # Sum everything
        for k,v in stats.items():
            if k not in nonbasekeys:
                stats['depth'] += len(v['baseq'])
                stats['bqualsum'] += sum(v['baseq'])
                if 'mapq' not in v:
                    stats['mqualsum'] += sum(v['baseq'])
                else:
                    stats['mqualsum'] += sum(v['mapq'])

        return stats

class StatsBase(Base):
    def setUp( self ):
        super( StatsBase, self ).setUp()
        # For reference
        #'G': { 'baseq': [40]*70 },
        #'A': { 'baseq': [40]*10 },
        #'C': { 'baseq': [40]*10 },
        #'T': { 'baseq': [40]*10 }
        #'depth': 100
        self.mock_stats()

class TestUnitBiasHQ(StatsBase):
    def _C( self, *args, **kwargs ):
        from base_caller import bias_hq
        return bias_hq( *args, **kwargs )

    def test_no_highquality( self ):
        r = self._C( self.stats )
        for k,v in self.stats.items():
            eq_( v, r[k] )

    def do_depth( self, stats, bias ):
        if stats is None:
            stats = self.stats
        # Depth should just increase by a factor of i
        r = self._C( stats, 1, bias )
        for k,v in stats.items():
            if k not in ('depth','mqualsum','bqualsum'):
                ebaseq = v['baseq']*int(bias)
                rbaseq = r[k]['baseq']
                eq_( ebaseq, r[k]['baseq'], 
                    "Len of base {} should be {} but got {}".format( k, len(ebaseq), len(rbaseq) )
                )
        # Verify depth is updated
        eq_( int(stats['depth']*bias), r['depth'] )

    def test_depth_is_updated( self ):
        # Verify that everything works when the bias is changed
        for i in (1,2,3.0):
            yield self.do_depth, None, i

    def test_biasth_works( self ):
        self.stats = {
            'A': {'baseq': [1,10,20,30,40,50,60] },
            'depth': 7
        }
        self.update_stats( self.stats )
        r = self._C( self.stats )
        eq_( self.stats['A']['baseq'] + [50,60], r['A']['baseq'] )
        r = self._C( self.stats, 15 )
        eq_( self.stats['A']['baseq'] + [20,30,40,50,60], r['A']['baseq'] )

    @raises(ValueError)
    def test_bias_is_zero( self ):
        # Non int and < 1 should raise error
        for i in (1, 0.9, 0, -1, 1.5):
            r = self._C( self.stats, 1, i )

class TestUnitMarkLQ(StatsBase):
    def _C( self, stats, minbq, mind ):
        from base_caller import mark_lq
        return mark_lq( stats, minbq, mind )

    def test_minq_lt( self ):
        self.stats['A']['baseq'] = [23,24,25,26]
        r = self._C( self.stats, 25, 1 )
        eq_( [23,24], r['?']['baseq'] )
        r = self._C( self.stats, 25, 100 )
        eq_( [23,24], r['N']['baseq'] )

    def test_removes_empty_bases( self ):
        self.stats['A']['baseq'] = [10]*10
        self.stats['C']['baseq'] = [10]*10
        r = self._C( self.stats, 25, 1 )
        assert 'A' not in r, 'A was not removed even though it had all < minq baseq'
        assert 'C' not in r, 'C was not removed even though it had all < minq baseq'
        eq_( [10]*20, r['?']['baseq'] )
        r = self._C( self.stats, 25, 100 )
        eq_( [10]*20, r['N']['baseq'] )

    def test_adds_n_single_base( self ):
        self.stats['A']['baseq'] = [10,10] + [30]*6 + [10,10]
        r = self._C( self.stats, 25, 1 )
        eq_( [10]*4, r['?']['baseq'] )
        r = self._C( self.stats, 25, 100 )
        eq_( [10]*4, r['N']['baseq'] )

    def test_no_n( self ):
        r = self._C( self.stats, 25, 1 )
        assert '?' not in r, '? was added to stats when it should not have'
        r = self._C( self.stats, 25, 100 )
        assert 'N' not in r, 'N was added to stats when it should not have'

class TestUnitCaller(Base):
    def setUp( self ):
        super( TestUnitCaller, self ).setUp()

    def _C( self, stats2, minbq, maxd, mind, minth ):
        from base_caller import caller
        return caller( stats2, minbq, maxd, mind, minth )

    def test_no_majority( self ):
        base_stats = {
            'A': { 'baseq': [25]*76 },
            'G': { 'baseq': [25]*12 },
            'C': { 'baseq': [25]*12 }
        }
        stats = self.make_stats( base_stats )
        # No > 0.8 majority, so call anything > 0.2 which should only be A
        eq_( ('A',76), self._C( stats, 25, 100000, 10, 0.8 ) )

    def test_calls_multiple_ambiguous( self ):
        # 33% A, 33% G, 33% C should end up V
        base_stats = {
            'A': { 'baseq': [25]*33 },
            'G': { 'baseq': [25]*33 },
            'C': { 'baseq': [25]*33 }
        }
        stats = self.make_stats( base_stats )
        # All bases are > 0.2 so return ambig AGC -> V
        eq_( ('V',99), self._C( stats, 25, 100000, 100, 0.8 ) )

    def test_calls_low_coverage_low_quality_n( self ):
        # 79% A, 21% Low Quality
        base_stats = {
            'A': { 'baseq': [25]*7 },
            '?': { 'baseq': [24]*3 }
        }
        stats = self.make_stats( base_stats )
        # Should trim off LQ G as they turn into ?
        eq_( ('N',3), self._C( stats, 25, 100000, 10, 0.8 ) )

    def test_calls_minth( self ):
        # 80% A, 20% G should be called A
        base_stats = {
            'A': { 'baseq': [40]*80 },
            'G': { 'baseq': [40]*20 }
        }
        stats = self.make_stats( base_stats )
        eq_( ('A',80), self._C( stats, 25, 100000, 10, 0.8 ) )

    def test_calls_specific_ambiguious( self ):
        # 79% A, 21% G should be called R
        base_stats = {
            'A': { 'baseq': [40]*79 },
            'G': { 'baseq': [40]*21 }
        }
        stats = self.make_stats( base_stats )
        eq_( ('R',100), self._C( stats, 25, 100000, 10, 0.8 ) )

    def test_example_1( self ):
        # Just testing an example from the test.vcf
        stats = {
            'bqualsum': 360.0,
			'mqualsum': 540.0,
			'depth': 9,
			'A': {
                'baseq': [40, 40, 40, 40, 40, 40, 40, 40],
                'mapq': [60, 60, 60, 60, 60, 60, 60, 60]
            },
            'C': {
                'baseq': [40],
                'mapq': [60]
            }
        }
        r = self._C( stats, 25, 100000, 10, 0.8 )
        eq_( ('A', 8), r )

    def test_example_2( self ):
        # Should call a D,3
        stats = {
			'bqualsum': 128.0,
			'mqualsum': 660.0,
			'depth': 3,
            'A': {
                'baseq': [40],
                'mapq': [60]
            },
			'G': {
                'baseq': [40],
                'mapq': [60]
            },
			'T': {
                'baseq': [40],
                'mapq': [60]
            }
        }
        r = self._C( stats, 25, 10000, 10, 0.8 )
        eq_( ('D',3), r )

class TestUnitCallOnPct(Base):
    def _C( self, stats, minth ):
        from base_caller import call_on_pct
        return call_on_pct( stats, minth )

    def test_called_ambigious( self ):
        base_stats = {
            'A': { 'baseq': [40]*80 },
            'C': { 'baseq': [40]*20 }
        }
        stats = self.make_stats( base_stats )
        r = self._C( stats, 0.9 )
        eq_( ('M',100), r )

    def test_gt_minth_gets_called( self ):
        base_stats = {
            'A': { 'baseq': [40]*80 },
            'C': { 'baseq': [40]*20 }
        }
        stats = self.make_stats( base_stats )
        r = self._C( stats, 0.8 )
        eq_( ('A',80), r )

    def test_calls_n( self ):
        base_stats = {
            'A': { 'baseq': [40]*10 },
            'C': { 'baseq': [40]*10 },
            'N': { 'baseq': [40]*80 }
        }
        stats = self.make_stats( base_stats )
        r = self._C( stats, 0.8 )
        eq_( ('N',80), r )

    def test_calls_correct_amb( self ):
        base_stats = {
            'A': { 'baseq': [40]*33 },
            'C': { 'baseq': [40]*33 },
            'G': { 'baseq': [40]*19 },
            'N': { 'baseq': [40]*15 }
        }
        stats = self.make_stats( base_stats )
        r = self._C( stats, 0.8 )
        eq_( ('M',66), r )

    def test_multiple_amb_called( self ):
        base_stats = {
            'A': { 'baseq': [40]*33 },
            'C': { 'baseq': [40]*33 },
            'G': { 'baseq': [40]*34 }
        }
        stats = self.make_stats( base_stats )
        r = self._C( stats, 0.8 )
        eq_( ('V',100), r )

    def test_no_majority( self ):
        base_stats = {
            'G': { 'baseq': [40]*76 },
            'A': { 'baseq': [40]*12 },
            'C': { 'baseq': [40]*12 }
        }
        stats = self.make_stats( base_stats )
        r = self._C( stats, 0.8 )
        eq_( ('G',76), r )

    def test_empty_stats_returns_n_zero( self ):
        stats = {}
        r = self._C( stats, 0.8 )
        eq_( ('N',0), r )

    def test_no_majority_returns_n_zero( self ):
        stats = {
            'A': { 'baseq': [40]*19 },
            'C': { 'baseq': [40]*19 },
            'G': { 'baseq': [40]*19 },
            'T': { 'baseq': [40]*19 },
            'N': { 'baseq': [40]*19 },
            '*': { 'baseq': [40]*1 },
            'depth': 100
        }
        r = self._C( stats, 0.8 )
        eq_( ('N', 100), r )

def eqs_( v1, v2, msg=None ):
    ''' Just run str on v1 and v2 and compare and use eq then '''
    if msg:
        eq_( str(v1), str(v2), msg )
    else:
        eq_( str(v1), str(v2) )

class TestUnitBlankVcfRows(Base):
    def _C( self, refname, refseq, curpos, lastpos, call='-' ):
        from base_caller import blank_vcf_rows
        return blank_vcf_rows( refname, refseq, curpos, lastpos, call )

    def test_gap_front( self ):
        # Should return 1 & 2
        r = self._C( 'ref', 'A'*10, 3, 0 )
        eq_( 1, r[0].POS )
        eq_( 2, r[1].POS )
    
    def test_gap_middle( self ):
        # Should return 4 & 5
        r = self._C( 'ref', 'A'*10, 6, 3 )
        eq_( 4, r[0].POS )
        eq_( 5, r[1].POS )
    
    def test_gap_end( self ): 
        # Same test case
        self.test_gap_middle()

    def test_single_gap( self ):
        # Should insert 4
        r = self._C( 'ref', 'A'*10, 5, 3 )
        eq_( 4, r[0].POS )

    def test_no_gap( self ):
        # Should return empty
        r = self._C( 'ref', 'A'*10, 3, 3 )
        r = self._C( 'ref', 'A'*10, 3, 4 )
        eq_( 0, len( r ) )

    def test_sets_call( self ):
        # Make sure that the call gets set correctly
        r = self._C( 'ref', 'A'*10, 6, 3, 'N' )
        eq_( 'N', r[0].INFO['CB'] )
        r = self._C( 'ref', 'A'*10, 6, 3, '-' )
        eq_( '-', r[0].INFO['CB'] )

    def test_sets_refbase( self ):
        # Should be C then G
        r = self._C( 'ref', 'ACGT', 4, 1, 'N' )
        eq_( 'C', r[0].REF )
        eq_( 'G', r[1].REF )
        # Gap at the end
        r = self._C( 'ref', 'ACGTA', 6, 4, 'N' )
        eq_( 'A', r[0].REF )

class TestUnitBlankVcfRow(Base):
    def _C( self, refname, refseq, pos, call ):
        from base_caller import blank_vcf_row
        return blank_vcf_row( refname, refseq, pos, call )

    def sets_( self, record, refbase, pos, call ):
        r = record
        eq_( pos, r.POS )
        eq_( refbase, r.REF )
        eq_( 'ref', r.CHROM )
        eq_( [], r.ALT )
        eq_( 0, r.INFO['DP'] )
        eq_( 0, r.INFO['RC'] )
        eq_( 0, r.INFO['RAQ'] )
        eq_( 0, r.INFO['PRC'] )
        eq_( call, r.INFO['CB'] )
        eq_( 0, r.INFO['CBD'] )

    def test_expected( self ):
        rseq = 'ACGT'*2
        for i in range(1, len(rseq)+1):
            r = self._C( 'ref', rseq, i, '-' )
            yield self.sets_, r, rseq[i-1], i, '-'

@patch('base_caller.MPileupColumn')
class TestUnitGenerateVcfRow(Base):
    def _C( self, mpilecol, refseq, minbq, maxd, mind, minth, biasth=50, bias=2 ):
        from base_caller import generate_vcf_row
        return generate_vcf_row( mpilecol, refseq, minbq, maxd, mind, minth, biasth, bias )

    def mock_stats( self ):
        base_stats = {
            'G': { 'baseq': [40]*70 },
            'A': { 'baseq': [40]*10 },
            'C': { 'baseq': [40]*10 },
            'T': { 'baseq': [40]*10 }
        }
        return self.make_stats( base_stats )

    def setup_mpileupcol( self, mpilemock, ref='ref', pos=1, stats=None ):
        ''' Sets up a mock mpileupcolumn object '''
        if stats is None:
            stats = self.mock_stats()
        mpilemock.ref = ref
        mpilemock.pos = pos
        mpilemock.base_stats.return_value = stats

    @timed(EXPECTED_TIME_PER_BASE)
    def test_runs_quickly( self, mpilecol ):
        self.setup_mpileupcol( mpilecol )
        r = self._C( mpilecol, 'ACGT'*100, 25, 1000, 10, 0.8 )

    def test_bias_works( self, mpilecol ):
        stats = self.mock_stats()
        # SHould keep depth at 100 and as it stands
        # ambiguous call which we can bias towards the N
        stats['G'] = {'baseq': [50]*60}
        stats['N'] = {'baseq': [40]*10}
        self.setup_mpileupcol( mpilecol, stats=stats )
        # biasth @ 50 will bias the G and bias @ 3 will make it
        # the majority at 180/220 = 82%
        r = self._C( mpilecol, 'ACGT'*100, 25, 1000, 10, 0.8, 50, 3 )
        eq_( 220, r.INFO['DP'] )
        eq_( 180, r.INFO['CBD'] )
        eq_( 'G', r.INFO['CB'] )

    def test_regionstr_not_1( self, mpilecol ):
        self.setup_mpileupcol( mpilecol )
        r = self._C( mpilecol, 'ACGT'*100, 25, 1000, 10, 0.8 )
        eqs_( 100, r.INFO['DP'] )

    def test_depth_set( self, mpilecol):
        self.setup_mpileupcol( mpilecol )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eqs_( 100, r.INFO['DP'] )

    def test_refstats_set( self, mpilecol):
        self.setup_mpileupcol( mpilecol )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eqs_( 10, r.INFO['RC'] )
        eqs_( 10, r.INFO['PRC'] )
        eqs_( 40, r.INFO['RAQ'] )

    def test_calls_amb( self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*9 },
            'C': { 'baseq': [38]*21 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        with patch( 'base_caller.info_stats' ) as info_stats:
            info_stats.return_value = {
                'AC':[60,21,9],
                'PAC':[60,21,10],
                'AAQ':[40,38,37],
                'bases': ['G','C','T']
            }
            r = self._C( mpilecol, 'A', 25, 1000, 100, 0.8 )
        eq_( [60,21,9], r.INFO['AC'] )
        eq_( [60,21,10], r.INFO['PAC'] )
        eq_( [40,38,37], r.INFO['AAQ'] )
        eq_( ['G','C','T'], r.ALT )
        eq_( 'S', r.INFO['CB'] )

    def test_alternatestats_same_order( self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*10 },
            'C': { 'baseq': [38]*20 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        with patch( 'base_caller.info_stats' ) as info_stats:
            info_stats.return_value = {
                'AC':[60,20,10],
                'PAC':[60,20,10],
                'AAQ':[40,38,37],
                'bases': ['G','C','T']
            }
            r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( [60,20,10], r.INFO['AC'] )
        eq_( [60,20,10], r.INFO['PAC'] )
        eq_( [40,38,37], r.INFO['AAQ'] )
        eq_( ['G','C','T'], r.ALT )

    def test_ensure_values_rounded( self, mpilecol):
        base_stats = {
            'G': { 'baseq': [38]*17+[39]*16 },
            'A': { 'baseq': [30]*33 },
            'T': { 'baseq': [33]*16 + [35]*17 } # Hopefully avg qual will be 33.33 and depth is 33% too
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        r = self._C( mpilecol, 'T', 25, 100, 10, 0.8 )
        eqs_( 33, r.INFO['RC'] )
        eqs_( 34, r.INFO['RAQ'] )
        eqs_( 33, r.INFO['PRC'] )
        bases = r.ALT
        ac = r.INFO['AC']
        aaq = r.INFO['AAQ']
        pac = r.INFO['PAC']
        altinfo = dict( zip( bases, zip( ac, aaq, pac ) ) )
        eqs_( 33, altinfo['G'][0] )
        eqs_( 38, altinfo['G'][1] )
        eqs_( 33, altinfo['G'][2] )
        eqs_( 33, altinfo['A'][0] )
        eqs_( 30, altinfo['A'][1] )
        eqs_( 33, altinfo['A'][2] )

    def test_single_alternate_base_stats( self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*90 },
            'A': { 'baseq': [40]*10 },
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( [90], r.INFO['AC'] )
        eq_( [90], r.INFO['PAC'] )
        eq_( [40], r.INFO['AAQ'] )

    def test_called_base( self, mpilecol):
        self.setup_mpileupcol( mpilecol )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( 'G', r.INFO['CB'] )
        eqs_( 70, r.INFO['CBD'] )

    def test_fields_set_multiple( self, mpilecol):
        self.setup_mpileupcol( mpilecol )
        # So we can check the ordering is correct
        info_stats = {'AC':[1,2,3],'AAQ':[4,5,6],'PAC':[7,8,9],'bases':['A','C','G']}
        with patch('base_caller.info_stats') as info_stats:
            info_stats.return_value = info_stats
            r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( 'ref', r.CHROM )
        eqs_( 1, r.POS )
        eqs_( 'A', r.REF )
        eq_( info_stats['bases'], r.ALT )

    def test_fields_set_single( self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*80 },
            'A': { 'baseq': [40]*10 }
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( 'ref', r.CHROM )
        eqs_( 1, r.POS )
        eq_( 'A', r.REF )
        eq_( ['G'], r.ALT )

    def test_refbase_not_in_stats( self, mpilecol):
        base_stats = {
            'N': { 'baseq': [40]*100 },
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( 'ref', r.CHROM )
        eqs_( 1, r.POS )
        eq_( 'A', r.REF )
        eq_( ['N'], r.ALT )
        eq_( 0, r.INFO['RC'] )
        eq_( 0, r.INFO['RAQ'] )
        eq_( 0, r.INFO['PRC'] )

    def test_depth_gte_mind( self, mpilecol):
        # Ensures that the A is turned into a ? and trimmed
        # as it is low quality
        stats = {
            'depth': 20,
            'A': { 'baseq': [1]*10 },
            'C': { 'baseq': [40]*10 },
            'mqualsum': 0,
            'bqualsum': 0
        }
        self.setup_mpileupcol( mpilecol, stats=stats )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( 'C', r.INFO['CB'] )
        eq_( 10, r.INFO['CBD'] )
        eq_( 0, r.INFO['RC'] )

    def test_depth_lt_mind( self, mpilecol):
        # Ensures that the A is turned into a ? and trimmed
        # as it is low quality
        stats = {
            'depth': 20,
            'A': { 'baseq': [1]*10 },
            'C': { 'baseq': [40]*10 },
            'mqualsum': 0,
            'bqualsum': 0
        }
        self.setup_mpileupcol( mpilecol, stats=stats )
        r = self._C( mpilecol, 'A', 25, 1000, 10, 0.8 )
        eq_( 'C', r.INFO['CB'] )
        eq_( 10, r.INFO['CBD'] )
        eq_( 0, r.INFO['RC'] )
        eq_( [10], r.INFO['AC'] )
        eq_( [100], r.INFO['PAC'] )

class TestUnitGenerateVCF(Base):
    # Hard to test each thing without generating sam files and vcf manually so
    # just going to let the integration tests do it...
    pass

class TestUnitInfoStats(Base):
    def setUp( self ):
        self.mock_stats2()

    def _C( self, stats2, refbase ):
        from base_caller import info_stats
        return info_stats( stats2, refbase )

    def mock_stats2( self ):
        ''' Every base and each base has length 20 with 40 quality so depth 100 '''
        s = { 'baseq': [40]*20 }
        self.stats2 = {'depth':0}
        for b in 'ACGNT':
            self.stats2['depth'] = len(s['baseq'])
            self.stats2[b] = s.copy()

    def test_empty_result( self ):
        self.stats2 = {
            'mqualsum': 4000,
            'bqualsum': 6000,
            'depth': 100,
            'A': {
                'baseq': [40]*100,
            }
        }
        r = self._C( self.stats2, 'A' )
        eq_( {'bases':[], 'AAQ': [], 'AC': [], 'PAC': []}, r )

    def test_excludes_keys( self ):
        self.stats2 = {}
        self.stats2['mqualsum'] = 1
        self.stats2['bqualsum'] = 1
        self.stats2['depth'] = 1
        self.stats2['A'] = 1
        r = self._C( self.stats2, 'A' )
        for k in ('AC','AAQ','PAC','bases'):
            eq_( [], r[k] )

    def test_ensure_expected( self ):
        ''' Order of the bases must be preserved '''
        self.stats2['A']['baseq'] = [40]*200
        self.stats2['C']['baseq'] = [39]*200
        self.stats2['G']['baseq'] = [38]*42
        self.stats2['N']['baseq'] = [37]*58
        self.stats2['T']['baseq'] = [36]*500
        self.stats2['depth'] = 1000
        r = self._C( self.stats2, 'A' )
        eq_( [200,42,58,500], r['AC'] )
        eq_( [20,4,6,50], r['PAC'] )
        eq_( [39,38,37,36], r['AAQ'] )
        eq_( ['C','G','N','T'], r['bases'] )

    def test_ensure_exclude_ref( self ):
        self.stats2['G']['baseq'] = [10]*20
        self.stats2['depth'] = 100
        r = self._C( self.stats2, 'G' )
        eq_( [40]*4, r['AAQ'] )
        eq_( ['A','C','N','T'], r['bases'] )

class BaseInty(Base):
    def print_files( self, f1, f2 ):
        print open(f1).read()
        print open(f2).read()

    def cmp_files( self, f1, f2 ):
        import subprocess
        try:
            assert exists( f1 ), "{} doesn't exist".format(f1)
            assert exists( f2 ), "{} doesn't exist".format(f2)
            subprocess.check_output( 'diff {} {}'.format(f1, f2), shell=True )
            return True
        except subprocess.CalledProcessError as e:
            print e
            print e.output
            print "Here is the generated file for comparison"
            print open(f2).read()
            return False

    def temp_bam( self, bam, bai ):
        ''' Copies given bam and bai to fixture temp dir '''
        import shutil
        tbam = join( self.tempdir, basename(self.bam) )
        tbai = tbam + '.bai'
        shutil.copy( self.bam, tbam )
        shutil.copy( self.bai, tbai )
        return tbam, tbai

class TestUnitGenerateVCF(BaseInty):
    def _C( self, bamfile, reffile, regionstr, vcf_output_file, minbq, maxd, mind=10, minth=10, biasth=50, bias=2, vcf_template=None ):
        from base_caller import generate_vcf, VCF_HEAD
        if vcf_template is None:
            template = VCF_HEAD.format(basename(bamfile))
        else:
            template = vcf_template
        return generate_vcf( bamfile, reffile, regionstr, vcf_output_file, 
                minbq, maxd, mind, minth, biasth, bias, vcf_template=template ) 

    def test_correct_vcf( self ):
        out_vcf = join( self.tempdir, 'out.vcf' )
        r = self._C( self.bam, self.ref, None, out_vcf, 25, 100, 10, 0.8 )
        assert self.cmp_files( self.vcf, out_vcf )
        eq_( out_vcf, r )

    def test_default_outfile( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        r = self._C( tbam, self.ref, None, None, 25, 100, 10, 0.8 )
        assert self.cmp_files( self.vcf, out_vcf )
        eq_( out_vcf, r )

    @timed(90)
    @attr('slow')
    def test_fullsample( self ):
        # The base directory of the fixture files
        fullsampledir = fsd = join( fixtures.THIS, 'fixtures', 'base_caller', 'fullsample' )
        # Input files
        ref = join( fsd, 'Den1__WestPac__1997.fasta' )
        bam = join( fsd, 'fullsample.bam' )
        bai = bam + '.bai'
        # Expected vcf
        vcf = join( fsd, 'fullsample.bam.vcf' )
        # Result vcf
        out_vcf = join( self.tempdir, basename(bam) + '.vcf' )
        print out_vcf
        assert out_vcf is not None
        # Run the base caller
        # Set bias to 10. Will probably be what will be used as default eventually
        out_vcf = self._C( bam, ref, None, out_vcf, 25, 100000, 10, 0.8, 50, 10 )
        # Compare the expected and result vcf
        assert self.cmp_files( vcf, out_vcf )


@attr('current')
class TestUnitMain(BaseInty):
    def _C( self, bamfile, reffile, vcf_output_file, regionstr=None, minbq=25, maxd=100000, mind=10, minth=0.8, biasth=50, bias=2 ):
        from base_caller import main
        args = Mock(
            bamfile=bamfile,
            reffile=reffile,
            vcf_output_file=vcf_output_file,
            regionstr=regionstr,
            minbq=minbq,
            maxd=maxd,
            mind=mind,
            minth=minth,
            biasth=biasth,
            bias=bias
        )        
        return main( args )

    def test_noregion_emptypileup( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        r = self._C( self.bam, self.ref, out_vcf, 'doesnotexist', 25, 100, 10, 0.8, 50, 2 )
        assert not self.cmp_files( self.vcf, out_vcf )

    def test_runs( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        r = self._C( self.bam, self.ref, out_vcf, '', 25, 100, 10, 0.8, 50, 2 )
        assert self.cmp_files( self.vcf, out_vcf )

class TestIntegrate(BaseInty):
    def _C( self, bamfile, reffile, vcf_output_file, regionstr=None, minbq=25, maxd=100000, mind=10, minth=0.8, biasth=50, bias=2 ):
        import subprocess
        script_path = join( dirname( dirname( __file__ ) ) )
        script_path = join( script_path, 'base_caller.py' )
        cmd = [script_path, bamfile, reffile]
        if regionstr:
            cmd += ['-r', regionstr]
        if vcf_output_file:
            cmd += ['-o', vcf_output_file]
        cmd += ['-minbq', minbq, '-maxd', maxd, '-mind', mind, '-minth', minth, '-biasth', biasth, '-bias', bias]
        cmd = [str(x) for x in cmd]
        #print cmd
        return subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE )

    @timed(EXPECTED_TOTAL_TIME)
    def test_exit_0_on_success( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        p = self._C( self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8 )
        eq_( 0, p.wait() )

    def test_outputs_correct_vcf_1( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        p = self._C( self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8 )
        o,e = p.communicate()
        assert self.cmp_files( self.vcf, out_vcf )

    def test_stdouterr_ok_and_filesmatch( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        p = self._C( self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8 )

        o,e = p.communicate()
        print "STDOUT"
        print o
        print "STDERR"
        print e
        assert e==o==''

        assert self.cmp_files( self.vcf, out_vcf )

    def test_nondefault_filesdiffer( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        p = self._C( self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8, 50, 10 )

        o,e = p.communicate()
        print "STDOUT"
        print o
        print "STDERR"
        print e
        assert e==o==''

        assert not self.cmp_files( self.vcf, out_vcf )
