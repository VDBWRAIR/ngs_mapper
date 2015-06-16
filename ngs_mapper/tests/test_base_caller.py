from imports import *
from ngs_mapper.samtools import InvalidRegionString

from ngs_mapper.base_caller import VCF_HEAD

# How long you expect each base position on the reference to take to process
EXPECTED_TIME_PER_BASE = 0.0015
# How long you expect for the entire fixtures/base_caller/test.bam to complete
EXPECTED_TOTAL_TIME = 0.6

class Base( common.BaseBaseCaller ):
    modulepath = 'ngs_mapper.base_caller'

    def mock_stats(self):
        s = { 'baseq': [40]*10, 'mapq': [60]*10 }
        x = {}
        for b in 'ATGC':
            x[b] = s.copy()
        self.stats = self.make_stats(x)

    def mock_mpileup_factory(self, **kwargs):
        '''
        Returns a function that can mock out ngs_mapper.samtools.mpileup
        
        pileupstart - start of the pileup that will be generated
        pileupend - end of the pileup that will be generated
        refdepth - depth of each reference base
        refbase - base that the reference will have for each pileup
        pileupbase - all the bases in the pileup will be this
        pileupqualstr - the character that will be used for all quals(mqual & bqual)
        '''
        refdepth = kwargs['refdepth']
        def get_mpileup_region(bamfile, regionstr, mind, minq, maxd):
            '''
            Mock ngs_mapper.samtools.mpileup region
            '''
            from ngs_mapper.samtools import parse_regionstring
            ref, rstart, rend = parse_regionstring(regionstr)
            # Have to pick min/max of these to restrict correctly
            pilestart = max(rstart, kwargs['pileupstart'])
            pileend = min(rend, kwargs['pileupend'])
            return [
                self._mock_pileup_str(
                    ref, i, kwargs['refbase'],
                    refdepth, kwargs['pileupbase']*refdepth,
                    kwargs['pileupqualstr']*refdepth, kwargs['pileupqualstr']*refdepth
               )
                for i in range(pilestart, pileend+1)
            ]
        return get_mpileup_region

    def make_stats(self, base_stats):
        stats = {}
        for k,v in base_stats.items():
            stats[k] = {}
            stats[k]['baseq'] = v['baseq']
            if 'mapq' not in stats[k]:
                stats[k]['mapq'] = v['baseq']
            else:
                stats[k]['mapq'] = v['mapq']

        self.update_stats(stats)

        return stats

    def update_stats(self, stats):
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

class Hpoly(Base):
    def make_hpoly_fasta(self, hlist):
        with open('ref.fasta', 'w') as fh:
            seq = ['X'] * hlist[-1][-1]
            for nucs, s, e in hlist:
                   seq[s-1:e] = nucs
            seq = ''.join(seq)
            fh.write('>ref\n{0}\n'.format(seq))
        return 'ref.fasta'

    def setUp(self):
        super(Hpoly, self).setUp()
        #AAATAAAATAAAAA
        self.hlist = [('AAA',1,3),('AAAA',5,8),('AAAAA',10,14)]
        self.ref = self.make_hpoly_fasta(self.hlist)
        self.seqs = SeqIO.index(self.ref, 'fasta')

class TestHpolyList(Hpoly):
    functionname = 'hpoly_list'

    def make_list(self, hpolys):
        # easier to compare
        r = {}
        for s in hpolys:
            r[s] = list(hpolys[s])
        return r

    def test_single_sequence(self):
        r = self._C(self.seqs, 3)
        r = self.make_list(r)
        eq_({'ref':self.hlist}, r)

    def test_minlength(self):
        r = self._C(self.seqs, 4)
        r = self.make_list(r)
        eq_({'ref':self.hlist[1:]}, r)

class TestIsHpoly(Hpoly):
    functionname = 'is_hpoly'

    def setUp( self ):
        super( TestIsHpoly, self ).setUp()
        from ngs_mapper.base_caller import hpoly_list
        self.hpoly = hpoly_list( self.seqs, 3 )

    def test_(self):
        ok_(self._C(self.hpoly, 'ref', 1))
        ok_(self._C(self.hpoly, 'ref', 2))
        ok_(self._C(self.hpoly, 'ref', 3))
        ok_(not self._C(self.hpoly, 'ref', 4))
        ok_(self._C(self.hpoly, 'ref', 5))
        ok_(self._C(self.hpoly, 'ref', 6))
        ok_(self._C(self.hpoly, 'ref', 7))
        ok_(self._C(self.hpoly, 'ref', 8))
        ok_(not self._C(self.hpoly, 'ref', 9))
        ok_(self._C(self.hpoly, 'ref', 10))
        ok_(self._C(self.hpoly, 'ref', 11))
        ok_(self._C(self.hpoly, 'ref', 12))
        ok_(self._C(self.hpoly, 'ref', 13))
        ok_(self._C(self.hpoly, 'ref', 14))

class StatsBase(Base):
    def setUp(self):
        super(StatsBase, self).setUp()
        # For reference
        #'G': { 'baseq': [40]*70 },
        #'A': { 'baseq': [40]*10 },
        #'C': { 'baseq': [40]*10 },
        #'T': { 'baseq': [40]*10 }
        #'depth': 100
        self.mock_stats()

class TestUnitBiasHQ(StatsBase):
    functionname = 'bias_hq'

    def test_no_highquality(self):
        r = self._C(self.stats)
        for k,v in self.stats.items():
            eq_(v, r[k])

    def do_depth(self, stats, bias):
        if stats is None:
            stats = self.stats
        # Depth should just increase by a factor of i
        r = self._C(stats, 1, bias)
        for k,v in stats.items():
            if k not in ('depth','mqualsum','bqualsum'):
                ebaseq = v['baseq']*int(bias)
                rbaseq = r[k]['baseq']
                eq_(ebaseq, r[k]['baseq'], 
                    "Len of base {0} should be {1} but got {2}".format(k, len(ebaseq), len(rbaseq))
               )
        # Verify depth is updated
        eq_(int(stats['depth']*bias), r['depth'])

    def test_depth_is_updated(self):
        # Verify that everything works when the bias is changed
        for i in (1,2,3.0):
            yield self.do_depth, None, i

    def test_biasth_works(self):
        self.stats = {
            'A': {'baseq': [1,10,20,30,40,50,60] },
            'depth': 7
        }
        self.update_stats(self.stats)
        r = self._C(self.stats, 50, 2)
        eq_(self.stats['A']['baseq'] + [50,60], r['A']['baseq'])
        r = self._C(self.stats, 15, 2)
        eq_(self.stats['A']['baseq'] + [20,30,40,50,60], r['A']['baseq'])

    @raises(ValueError)
    def test_bias_is_zero(self):
        # Non int and < 1 should raise error
        for i in (1, 0.9, 0, -1, 1.5):
            r = self._C(self.stats, 1, i)

class TestUnitMarkLQ(StatsBase):
    functionname = 'mark_lq'

    def test_bias_reference(self):
        # If less than mind then don't call N if base is == ref base
        # ideally the following would become T
        # Issue #395
        base_stats = {
            'T': {
                'baseq': [16, 17, 17, 19, 21, 29, 32, 35, 37]
            }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 25, 10, 'T')
        eq_(9, len(r['T']['baseq']))

        # Now add an A as well just to make sure. Also puts us >= mind of 10
        base_stats['A'] = {'baseq':[16]}
        stats = self.make_stats(base_stats)
        r = self._C(stats, 25, 10, 'T')
        eq_(4, len(r['T']['baseq']))
        eq_(6, len(r['?']['baseq']))

    def test_minq_lt(self):
        self.stats['A']['baseq'] = [23,24,25,26]
        r = self._C(self.stats, 25, 1, 'G')
        eq_([23,24], r['?']['baseq'])
        r = self._C(self.stats, 25, 100, 'G')
        eq_([23,24], r['N']['baseq'])

    def test_removes_empty_bases(self):
        self.stats['A']['baseq'] = [10]*10
        self.stats['C']['baseq'] = [10]*10
        r = self._C(self.stats, 25, 1, 'G')
        assert 'A' not in r, 'A was not removed even though it had all < minq baseq'
        assert 'C' not in r, 'C was not removed even though it had all < minq baseq'
        eq_([10]*20, r['?']['baseq'])
        r = self._C(self.stats, 25, 100, 'G')
        eq_([10]*20, r['N']['baseq'])

    def test_adds_n_single_base(self):
        # A - Depth 10, AQ - 20
        self.stats['A']['baseq'] = [10,10] + [30]*6 + [10,10]
        r = self._C(self.stats, 25, 1, 'G')
        eq_([10]*4, r['?']['baseq'])
        r = self._C(self.stats, 25, 100, 'G')
        eq_([10]*4, r['N']['baseq'])

    def test_ensure_depth_threshold(self):
        # Make sure if the depth is == mind everything still works
        # Should trigger all N
        stats = self.make_stats({
            'A': {'baseq':[20]*10}
        })
        # A's should all get turned to ? because mind == depth is high coverage
        r = self._C(stats, 25, 10, 'G')
        eq_(10, len(r['?']['baseq']))

    def test_no_n(self):
        r = self._C(self.stats, 25, 1, 'G')
        assert '?' not in r, '? was added to stats when it should not have'
        r = self._C(self.stats, 25, 100, 'G')
        assert 'N' not in r, 'N was added to stats when it should not have'

class TestUnitCaller(Base):
    functionname = 'caller'

    def setUp(self):
        super(TestUnitCaller, self).setUp()

    def test_no_majority(self):
        # No > 0.8 majority, so call anything > 0.2 which should only be A
        base_stats = {
            'A': { 'baseq': [25]*76 },
            'G': { 'baseq': [25]*12 },
            'C': { 'baseq': [25]*12 }
        }
        stats = self.make_stats(base_stats)
        eq_(('A',76), self._C(stats, 25, 100000, 10, 0.8))

    def test_calls_multiple_ambiguous(self):
        # 33% A, 33% G, 33% C should end up V
        base_stats = {
            'A': { 'baseq': [25]*33 },
            'G': { 'baseq': [25]*33 },
            'C': { 'baseq': [25]*33 }
        }
        stats = self.make_stats(base_stats)
        # All bases are > 0.2 so return ambig AGC -> V
        eq_(('V',99), self._C(stats, 25, 100000, 100, 0.8))

    def test_calls_low_coverage_low_quality_n(self):
        # 79% A, 21% Low Quality
        base_stats = {
            'A': { 'baseq': [25]*7 },
            '?': { 'baseq': [24]*3 }
        }
        stats = self.make_stats(base_stats)
        # Should trim off LQ G as they turn into ?
        eq_(('N',3), self._C(stats, 25, 100000, 10, 0.8))

    def test_calls_minth(self):
        # 80% A, 20% G should be called A
        base_stats = {
            'A': { 'baseq': [40]*80 },
            'G': { 'baseq': [40]*20 }
        }
        stats = self.make_stats(base_stats)
        eq_(('A',80), self._C(stats, 25, 100000, 10, 0.8))

    def test_calls_specific_ambiguious(self):
        # 79% A, 21% G should be called R
        base_stats = {
            'A': { 'baseq': [40]*79 },
            'G': { 'baseq': [40]*21 }
        }
        stats = self.make_stats(base_stats)
        eq_(('R',100), self._C(stats, 25, 100000, 10, 0.8))

    def test_example_1(self):
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
        r = self._C(stats, 25, 100000, 10, 0.8)
        eq_(('A', 8), r)

    def test_example_2(self):
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
        r = self._C(stats, 25, 10000, 10, 0.8)
        eq_(('D',3), r)

class TestUnitCallOnPct(Base):
    functionname = 'call_on_pct'

    def test_called_ambigious(self):
        base_stats = {
            'A': { 'baseq': [40]*80 },
            'C': { 'baseq': [40]*20 }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 0.9)
        eq_(('M',100), r)

    def test_gt_minth_gets_called(self):
        base_stats = {
            'A': { 'baseq': [40]*80 },
            'C': { 'baseq': [40]*20 }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 0.8)
        eq_(('A',80), r)

    def test_calls_n(self):
        base_stats = {
            'A': { 'baseq': [40]*10 },
            'C': { 'baseq': [40]*10 },
            'N': { 'baseq': [40]*80 }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 0.8)
        eq_(('N',80), r)

    def test_calls_correct_amb(self):
        base_stats = {
            'A': { 'baseq': [40]*33 },
            'C': { 'baseq': [40]*33 },
            'G': { 'baseq': [40]*19 },
            'N': { 'baseq': [40]*15 }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 0.8)
        eq_(('M',66), r)

    def test_multiple_amb_called(self):
        base_stats = {
            'A': { 'baseq': [40]*33 },
            'C': { 'baseq': [40]*33 },
            'G': { 'baseq': [40]*34 }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 0.8)
        eq_(('V',100), r)

    def test_N_and_other_base_gt_threshold(self):
        base_stats = {
            'A': { 'baseq': [40]*42 },
            'N': { 'baseq': [40]*42 },
            'G': { 'baseq': [40]*16 }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 0.8)
        eq_(('N',84), r)

    def test_no_majority(self):
        base_stats = {
            'G': { 'baseq': [40]*76 },
            'A': { 'baseq': [40]*12 },
            'C': { 'baseq': [40]*12 }
        }
        stats = self.make_stats(base_stats)
        r = self._C(stats, 0.8)
        eq_(('G',76), r)

    def test_empty_stats_returns_gap_zero(self):
        stats = {}
        r = self._C(stats, 0.8)
        eq_(('-',0), r)

    def test_no_majority_returns_n_zero(self):
        stats = {
            'A': { 'baseq': [40]*19 },
            'C': { 'baseq': [40]*19 },
            'G': { 'baseq': [40]*19 },
            'T': { 'baseq': [40]*19 },
            'N': { 'baseq': [40]*19 },
            '*': { 'baseq': [40]*1 },
            'depth': 100
        }
        r = self._C(stats, 0.8)
        eq_(('N', 100), r)

def eqs_(v1, v2, msg=None):
    ''' Just run str on v1 and v2 and compare and use eq then '''
    if msg:
        eq_(str(v1), str(v2), msg)
    else:
        eq_(str(v1), str(v2))

class TestUnitBlankVcfRows(Base):
    functionname = 'blank_vcf_rows'

    def test_gap_front(self):
        # Should return 1 & 2
        r = self._C('ref', 'A'*10, 0, 3)
        eq_(1, r[0].POS)
        eq_(2, r[1].POS)
    
    def test_gap_middle(self):
        # Should return 4 & 5
        r = self._C('ref', 'A'*10, 3, 6)
        eq_(4, r[0].POS)
        eq_(5, r[1].POS)
    
    def test_gap_end(self): 
        r = self._C('ref', 'A'*10, 9, 11)
        eq_(10, r[0].POS)

    def test_single_gap(self):
        # Should insert 4
        r = self._C('ref', 'A'*10, 3, 5)
        eq_(4, r[0].POS)

    def test_no_gap(self):
        # Should return empty
        r = self._C('ref', 'A'*10, 3, 3)
        r = self._C('ref', 'A'*10, 4, 3)
        eq_(0, len(r))

    def test_sets_call(self):
        # Make sure that the call gets set correctly
        r = self._C('ref', 'A'*10, 3, 6, 'N')
        eq_('N', r[0].INFO['CB'])
        r = self._C('ref', 'A'*10, 3, 6, '-')
        eq_('-', r[0].INFO['CB'])

    def test_sets_refbase(self):
        # Should be C then G
        r = self._C('ref', 'ACGT', 1, 4, 'N')
        eq_('C', r[0].REF)
        eq_('G', r[1].REF)
        # Gap at the end
        r = self._C('ref', 'ACGTA', 4, 6, 'N')
        eq_('A', r[0].REF)

class TestUnitBlankVcfRow(Base):
    functionname = 'blank_vcf_row'

    def sets_(self, record, refbase, pos, call):
        r = record
        eq_(pos, r.POS)
        eq_(refbase, r.REF)
        eq_('ref', r.CHROM)
        eq_('.', r.ALT)
        eq_(0, r.INFO['DP'])
        eq_(0, r.INFO['RC'])
        eq_(0, r.INFO['RAQ'])
        eq_(0, r.INFO['PRC'])
        eq_(call, r.INFO['CB'])
        eq_(0, r.INFO['CBD'])

    def test_expected(self):
        rseq = 'ACGT'*2
        for i in range(1, len(rseq)+1):
            r = self._C('ref', rseq, i, '-')
            yield self.sets_, r, rseq[i-1], i, '-'

class MpileBase(Base):
    def mock_stats(self):
        base_stats = {
            'G': { 'baseq': [40]*70 },
            'A': { 'baseq': [40]*10 },
            'C': { 'baseq': [40]*10 },
            'T': { 'baseq': [40]*10 }
        }
        return self.make_stats(base_stats)

    def setup_mpileupcol(self, mpilemock, ref='ref', pos=1, stats=None):
        ''' Sets up a mock mpileupcolumn object '''
        if stats is None:
            stats = self.mock_stats()
        mpilemock.ref = ref
        mpilemock.pos = pos
        mpilemock.base_stats.return_value = stats


@patch('ngs_mapper.base_caller.MPileupColumn')
class TestUnitPileStats(MpileBase):
    functionname = 'pile_stats'

    def test_bias_ref_and_bias_hq(self, mpilecol):
        # Should bias 2 reads and bias ref
        base_stats = {
            'T': {'baseq': [24, 25, 46, 56, 58, 38]}
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'T', 25, 10, 50, 10)
        eq_(23, r['depth'])

@patch('ngs_mapper.base_caller.MPileupColumn')
class TestUnitGenerateVcfRow(MpileBase):
    functionname = 'generate_vcf_row'

    def test_issue_1012_calls_alt_or_ref(self, mpilecol):
        base_stats = {
            'C': {'baseq':[40]*50+[1]*10},
            '*': {'baseq':[40]*45},
            'G': {'baseq':[40]*1},
            'A': {'baseq':[40]*1},
            'N': {'baseq':[0]*1},
            'T': {'baseq':[20]*1}
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
    
        # Ensure Reference is called
        r = self._C(mpilecol, 'C', 25, 1000, 10, 0.8, 50, 10)
        eq_('C', r.INFO['CB'])
        eq_(50, r.INFO['CBD'])

        # Ensure Alt is called
        r = self._C(mpilecol, 'G', 25, 1000, 10, 0.8, 50, 10)
        eq_('C', r.INFO['CB'])
        eq_(50, r.INFO['CBD'])

    def test_issue_1012_calls_n(self, mpilecol):
        base_stats = {
            'C': {'baseq':[40]*45},
            '*': {'baseq':[40]*50},
            'G': {'baseq':[40]*1},
            'A': {'baseq':[40]*1},
            'N': {'baseq':[0]*1},
            'T': {'baseq':[20]*1}
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'C', 25, 1000, 10, 0.8, 50, 10)
        eq_('N', r.INFO['CB'])
        eq_(50, r.INFO['CBD'])

    def test_issue_1012_alt_count_sum_50_ref_50(self, mpilecol):
        # Gap and other base depth == 50% as well, but still should call C
        base_stats = {
            'C': {'baseq':[40]*50},
            '*': {'baseq':[40]*25},
            'G': {'baseq':[40]*25},
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'C', 25, 1000, 10, 0.8, 50, 10)
        eq_('C', r.INFO['CB'])
        eq_(50, r.INFO['CBD'])

    def test_thing(self, mpilecol):
        base_stats = {
            'T': {'baseq':[29, 35, 37, 32]}
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'T', 25, 1000, 10, 0.8, 50, 10)
        eq_('T', r.INFO['CB'])

    def test_low_coverage_low_quality_ref(self, mpilecol):
        base_stats = {
            'T': {
                'baseq': [16, 17, 17, 19, 21, 29, 32, 35, 37]
            }
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'T', 25, 1000, 10, 0.8, 50, 10)
        eq_('T', r.INFO['CB'])

    def test_highcov_lowqual(self, mpilecol):
        # Greater than mind but all lq bases should turn into a N
        self.setup_mpileupcol(mpilecol)
        # Set the low quality threshold to 41 to trigger all bases as low quality
        r = self._C(mpilecol, 'G', 41, 1000, 10, 0.8)
        eq_('N', r.INFO['CB'])
        # All values should be 0
        for k in ('DP','RC','RAQ','PRC','CBD'):
            eq_(0, r.INFO[k])

    def test_reference_lowercase_dna(self, mpilecol):
        # Issue #145
        # Reference bases are lowercase should be converted to uppercase
        self.setup_mpileupcol(mpilecol)
        r = self._C(mpilecol, 'g', 25, 1000, 10, 0.8)
        eq_('G', r.REF)
        eq_(70, r.INFO['RC'])
        eq_(70, r.INFO['PRC'])
        eq_(40, r.INFO['RAQ'])

    def test_no_alternates(self, mpilecol):
        stats = self.mock_stats()
        del stats['A']
        del stats['C']
        del stats['T']
        stats['depth'] = 70
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'G', 25, 1000, 10, 0.8)
        eq_('.', r.ALT)

    @timed(EXPECTED_TIME_PER_BASE)
    def test_runs_quickly(self, mpilecol):
        self.setup_mpileupcol(mpilecol)
        r = self._C(mpilecol, 'ACGT'*100, 25, 1000, 10, 0.8)

    def test_bias_works(self, mpilecol):
        stats = self.mock_stats()
        # SHould keep depth at 100 and as it stands
        # ambiguous call which we can bias towards the N
        stats['G'] = {'baseq': [50]*60}
        stats['N'] = {'baseq': [40]*10}
        self.setup_mpileupcol(mpilecol, stats=stats)
        # biasth @ 50 will bias the G and bias @ 3 will make it
        # the majority at 180/220 = 82%
        r = self._C(mpilecol, 'ACGT'*100, 25, 1000, 10, 0.8, 50, 3)
        eq_(220, r.INFO['DP'])
        eq_(180, r.INFO['CBD'])
        eq_('G', r.INFO['CB'])

    def test_regionstr_not_1(self, mpilecol):
        self.setup_mpileupcol(mpilecol)
        r = self._C(mpilecol, 'ACGT'*100, 25, 1000, 10, 0.8)
        eqs_(100, r.INFO['DP'])

    def test_depth_set(self, mpilecol):
        self.setup_mpileupcol(mpilecol)
        r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eqs_(100, r.INFO['DP'])

    def test_refstats_set(self, mpilecol):
        self.setup_mpileupcol(mpilecol)
        r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eqs_(10, r.INFO['RC'])
        eqs_(10, r.INFO['PRC'])
        eqs_(40, r.INFO['RAQ'])

    def test_calls_amb(self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*9 },
            'C': { 'baseq': [38]*21 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        with patch('ngs_mapper.base_caller.info_stats' ) as info_stats:
            info_stats.return_value = {
                'AC':[60,21,9],
                'PAC':[60,21,10],
                'AAQ':[40,38,37],
                'bases': ['G','C','T']
            }
            r = self._C(mpilecol, 'A', 25, 1000, 100, 0.8)
        eq_([60,21,9], r.INFO['AC'])
        eq_([60,21,10], r.INFO['PAC'])
        eq_([40,38,37], r.INFO['AAQ'])
        eq_(['G','C','T'], r.ALT)
        eq_('S', r.INFO['CB'])

    def test_alternatestats_same_order(self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*10 },
            'C': { 'baseq': [38]*20 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        self.setup_mpileupcol( mpilecol, stats=stats )
        with patch('ngs_mapper.base_caller.info_stats' ) as info_stats:
            info_stats.return_value = {
                'AC':[60,20,10],
                'PAC':[60,20,10],
                'AAQ':[40,38,37],
                'bases': ['G','C','T']
            }
            r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eq_([60,20,10], r.INFO['AC'])
        eq_([60,20,10], r.INFO['PAC'])
        eq_([40,38,37], r.INFO['AAQ'])
        eq_(['G','C','T'], r.ALT)

    def test_ensure_values_rounded(self, mpilecol):
        base_stats = {
            'G': { 'baseq': [38]*17+[39]*16 },
            'A': { 'baseq': [30]*33 },
            'T': { 'baseq': [33]*16 + [35]*17 } # Hopefully avg qual will be 33.33 and depth is 33% too
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'T', 25, 100, 10, 0.8)
        eqs_(33, r.INFO['RC'])
        eqs_(34, r.INFO['RAQ'])
        eqs_(33, r.INFO['PRC'])
        bases = r.ALT
        ac = r.INFO['AC']
        aaq = r.INFO['AAQ']
        pac = r.INFO['PAC']
        altinfo = dict(zip(bases, zip(ac, aaq, pac)))
        eqs_(33, altinfo['G'][0])
        eqs_(38, altinfo['G'][1])
        eqs_(33, altinfo['G'][2])
        eqs_(33, altinfo['A'][0])
        eqs_(30, altinfo['A'][1])
        eqs_(33, altinfo['A'][2])

    def test_single_alternate_base_stats(self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*90 },
            'A': { 'baseq': [40]*10 },
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eq_([90], r.INFO['AC'])
        eq_([90], r.INFO['PAC'])
        eq_([40], r.INFO['AAQ'])

    def test_called_base(self, mpilecol):
        self.setup_mpileupcol(mpilecol)
        r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eq_('G', r.INFO['CB'])
        eqs_(70, r.INFO['CBD'])

    def test_fields_set_multiple(self, mpilecol):
        self.setup_mpileupcol(mpilecol)
        # So we can check the ordering is correct
        info_stats = {'AC':[1,2,3],'AAQ':[4,5,6],'PAC':[7,8,9],'bases':['A','C','G']}
        with patch('ngs_mapper.base_caller.info_stats') as info_stats:
            info_stats.return_value = info_stats
            r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eq_('ref', r.CHROM)
        eqs_(1, r.POS)
        eqs_('A', r.REF)
        eq_(info_stats['bases'], r.ALT)

    def test_fields_set_single(self, mpilecol):
        base_stats = {
            'G': { 'baseq': [40]*80 },
            'A': { 'baseq': [40]*10 }
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eq_('ref', r.CHROM)
        eqs_(1, r.POS)
        eq_('A', r.REF)
        eq_(['G'], r.ALT)

    def test_refbase_not_in_stats(self, mpilecol):
        base_stats = {
            'N': { 'baseq': [40]*100 },
        }
        stats = self.make_stats(base_stats)
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eq_('ref', r.CHROM)
        eqs_(1, r.POS)
        eq_('A', r.REF)
        eq_(['N'], r.ALT)
        eq_(0, r.INFO['RC'])
        eq_(0, r.INFO['RAQ'])
        eq_(0, r.INFO['PRC'])

    def test_depth_gte_mind(self, mpilecol):
        # Ensures that the A is turned into a ? and trimmed
        # as it is low quality
        stats = {
            'depth': 20,
            'A': { 'baseq': [1]*10 },
            'C': { 'baseq': [40]*10 },
            'mqualsum': 0,
            'bqualsum': 0
        }
        self.setup_mpileupcol(mpilecol, stats=stats)
        r = self._C(mpilecol, 'A', 25, 1000, 10, 0.8)
        eq_('C', r.INFO['CB'])
        eq_(10, r.INFO['CBD'])
        eq_(0, r.INFO['RC'])

    def test_depth_lt_mind(self, mpilecol):
        # Ensures that the A is turned into a ? and trimmed
        # as it is low quality
        stats = {
            'depth': 20,
            'A': { 'baseq': [1]*10 },
            'C': { 'baseq': [40]*10 },
            'mqualsum': 0,
            'bqualsum': 0
        }
        self.setup_mpileupcol(mpilecol, stats=stats)
        # Should ensure that >= mind is working since 20 == stats['depth']
        r = self._C(mpilecol, 'A', 25, 1000, 20, 0.8)
        eq_('C', r.INFO['CB'])
        eq_(10, r.INFO['CBD'])
        eq_(0, r.INFO['RC'])
        eq_([10], r.INFO['AC'])
        eq_([100], r.INFO['PAC'])
        eq_(['C'], r.ALT)

class TestUnitGenerateVCF(Base):
    # Hard to test each thing without generating sam files and vcf manually so
    # just going to let the integration tests do it...
    pass

class TestUnitInfoStats(Base):
    functionname = 'info_stats'

    def setUp(self):
        self.mock_stats2()

    def mock_stats2(self):
        ''' Every base and each base has length 20 with 40 quality so depth 100 '''
        s = { 'baseq': [40]*20 }
        self.stats2 = {'depth':0}
        for b in 'ACGNT':
            self.stats2['depth'] = len(s['baseq'])
            self.stats2[b] = s.copy()

    def test_empty_result(self):
        self.stats2 = {
            'mqualsum': 4000,
            'bqualsum': 6000,
            'depth': 100,
            'A': {
                'baseq': [40]*100,
            }
        }
        r = self._C(self.stats2, 'A')
        eq_({'bases':[], 'AAQ': [], 'AC': [], 'PAC': []}, r)

    def test_excludes_keys(self):
        self.stats2 = {}
        self.stats2['mqualsum'] = 1
        self.stats2['bqualsum'] = 1
        self.stats2['depth'] = 1
        self.stats2['A'] = 1
        r = self._C(self.stats2, 'A')
        for k in ('AC','AAQ','PAC','bases'):
            eq_([], r[k])

    def test_ensure_expected(self):
        ''' Order of the bases must be preserved '''
        self.stats2['A']['baseq'] = [40]*200
        self.stats2['C']['baseq'] = [39]*200
        self.stats2['G']['baseq'] = [38]*42
        self.stats2['N']['baseq'] = [37]*58
        self.stats2['T']['baseq'] = [36]*500
        self.stats2['depth'] = 1000
        r = self._C(self.stats2, 'A')
        eq_([200,42,58,500], r['AC'])
        eq_([20,4,6,50], r['PAC'])
        eq_([39,38,37,36], r['AAQ'])
        eq_(['C','G','N','T'], r['bases'])

    def test_ensure_exclude_ref(self):
        self.stats2['G']['baseq'] = [10]*20
        self.stats2['depth'] = 100
        r = self._C(self.stats2, 'G')
        eq_([40]*4, r['AAQ'])
        eq_(['A','C','N','T'], sorted(r['bases']))

class BaseInty(Base):
    def print_files(self, f1, f2):
        print open(f1).read()
        print open(f2).read()

    def cmp_files(self, f1, f2):
        import subprocess
        try:
            assert exists(f1), "{0} doesn't exist".format(f1)
            assert exists(f2), "{0} doesn't exist".format(f2)
            print subprocess.check_output('diff {0} {1}'.format(f1, f2), shell=True)
            return True
        except subprocess.CalledProcessError as e:
            print e
            print e.output
            print "Here is the generated file for comparison"
            print open(f2).read()
            return False

    def _iter_two_vcf(self, vcf1, vcf2, startpos=1, endpos=None):
        '''
        Skip the header and iterate same position in each file
        Assume both vcf have same positions in them and same references

        Can specify startpos and endpos to compare
        '''
        if isinstance(vcf1,str):
            fh1 = open(vcf1)
        else:
            fh1 = vcf1
        if isinstance(vcf2,str):
            fh2 = open(vcf2)
        else:
            fh2 = vcf2

        # Burn off headers
        for line in fh1:
            if line.startswith('#'):
                continue
            break
        for line in fh2:
            if line.startswith('#'):
                continue
            break

        from itertools import izip
        for v1, v2 in izip(fh1, fh2):
            yield (
                self._split_vcfrow(v1),
                self._split_vcfrow(v2)
           )

    def _split_vcfrow(self, row):
        ref,pos,id,refbase,altbase,qual,filter,info = row.split('\t')
        info = re.findall('(\w+)=([-0-9a-zA-z_,]+)', info)
        info = dict(info)
        return {
            'CHROM':ref,
            'POS':pos, 
            'ID': id,
            'REF': refbase,
            'ALT': altbase,
            'QUAL': qual,
            'FILTER': filter,
            'INFO': info
        }

    def temp_bam(self, bam, bai):
        ''' Copies given bam and bai to fixture temp dir '''
        import shutil
        tbam = join(self.tempdir, basename(self.bam))
        tbai = tbam + '.bai'
        shutil.copy(self.bam, tbam)
        shutil.copy(self.bai, tbai)
        return tbam, tbai

class TestGenerateVCF(BaseInty):
    functionname = 'generate_vcf'

    def setUp(self):
        super(TestGenerateVCF,self).setUp()

        from ngs_mapper.base_caller import VCF_HEAD
        self.vcf_head = VCF_HEAD.format(basename('test.bam'))

    @raises(InvalidRegionString)
    @patch('ngs_mapper.base_caller.SeqIO')
    @patch('ngs_mapper.base_caller.mpileup')
    def test_raises_exception_regionstring_invalid(self, *args):
        self._C('test.bam', 'test.ref', None, 'out.vcf', 0, 0, 10, 0.8, 50, 10, VCF_HEAD, False)

    @patch('ngs_mapper.base_caller.SeqIO')
    @patch('ngs_mapper.base_caller.mpileup')
    def test_regionstr_lt0_and_gt_reflen(self, mmpileup, mseqio):
        mseqio.index.return_value = {'Ref1':MagicMock(seq='A'*10,id='Ref1')}
        mmpileup.side_effect = self.mock_mpileup_factory(
            pileupstart=1,
            pileupend=10,
            refdepth=10,
            refbase='A',
            pileupbase='A',
            pileupqualstr='I'
        )
        r = self._C('test.bam', 'test.ref', 'Ref1:0-30', 'out.vcf', 0, 0, 10, 0.8, 50, 10, VCF_HEAD, False)

    @patch('ngs_mapper.base_caller.SeqIO')
    @patch('ngs_mapper.base_caller.mpileup')
    def test_ref_in_bam_only_contains_some_bases(self, mmpileup, mseqio):
        reflen = 8
        refdepth = 10
        # Pileup only returns pos 5, 6, and 7
        # The ends need to be filled by generate_vcf
        mmpileup.side_effect = self.mock_mpileup_factory(
            pileupstart=5,
            pileupend=7,
            refdepth=refdepth,
            refbase='A',
            pileupbase='A',
            pileupqualstr='I'
       )
        refrec = Mock(seq='A'*reflen,id='Ref1')
        mseqio.parse.return_value = refrec
        mseqio.index.return_value = {'Ref1':refrec}


        def countvcf(start, expectedlinecount):
            linecount = 0
            i = start
            for line in open('out.vcf'):
                if not line.startswith('#'):
                    linecount += 1
                    print line
                    _line = line.split('\t')
                    eq_(str(i), _line[1])
                    if i in (1,2,3,4,8):
                        eq_('DP=0', _line[-1].split(';')[0])
                    else:
                        eq_('DP=10', _line[-1].split(';')[0])
                    i += 1
            eq_(expectedlinecount, linecount)

        # regionstr set so all bases should be returned
        r = self._C('test.bam','test.ref', 'Ref1:1-4', 'out.vcf', 25, 100, 10, 0.8)
        countvcf(1, 4)
        r = self._C('test.bam','test.ref', 'Ref1:5-8', 'out.vcf', 25, 100, 10, 0.8)
        countvcf(5, 4)
        r = self._C('test.bam','test.ref', 'Ref1:5-7', 'out.vcf', 25, 100, 10, 0.8)
        countvcf(5,3)

    @timed(90)
    @attr('slow')
    def test_fullsample_correct_called_bases_hpoly(self):
        # The base directory of the fixture files
        fullsampledir = fsd = join(fixtures.THIS, 'fixtures', 'base_caller', 'fullsample')
        # Input files
        ref = join(fsd, 'Den1__WestPac__1997.fasta')
        bam = join(fsd, 'fullsample.bam')
        bai = bam + '.bai'
        # Expected vcf
        vcf = join(fsd, 'fullsample.bam.vcf')
        # Result vcf
        out_vcf = join(self.tempdir, basename(bam) + '.vcf')
        print out_vcf
        assert out_vcf is not None
        # Run the base caller
        # Set bias to 10. Will probably be what will be used as default eventually
        out_vcf = self._C(bam, ref, 'Den1/U88535_1/WestPac/1997/Den1_1:1-15000', out_vcf, 25, 100000, 10, 0.8, 50, 10)

        # Compare the expected and result vcf
        for evcf, rvcf in self._iter_two_vcf(vcf, out_vcf):
            ecb = evcf['INFO']['CB']
            rcb = rvcf['INFO']['CB']
            eq_(evcf['REF'],rvcf['REF'])
            eq_(evcf['POS'],rvcf['POS'])
            print evcf['POS']
            eq_(ecb, rcb)

class TestGenerateVcfMultithreaded(BaseInty):
    functionname = 'generate_vcf_multithreaded'

    @patch('ngs_mapper.base_caller.os')
    @patch('ngs_mapper.base_caller.SeqIO')
    @patch('ngs_mapper.base_caller.mpileup')
    def test_correct_amount_lines(self, mmpileup, mseqio, mos):
        reflen = 10
        numrefs = 3
        refdepth = 100
        threads = 4
        # Each reference is 1 million bases long
        ref1 = Mock(seq='A'*reflen,id='Ref1')
        ref2 = Mock(seq='T'*reflen,id='Ref2')
        ref3 = Mock(seq='G'*reflen,id='Ref3')
        mseqio.index.return_value = {'Ref1':ref1, 'Ref2': ref2, 'Ref3': ref3}
        mseqio.parse.return_value = iter([ref1, ref2, ref3])

        mmpileup.side_effect = self.mock_mpileup_factory(
            pileupstart=1,
            pileupend=reflen,
            refdepth=refdepth,
            refbase='A',
            pileupbase='A',
            pileupqualstr='I'
       )

        out_vcf = self._C('in.bam', 'in.ref', 'out.vcf', 25, 100000, 10, 0.8, 50, 10, threads)

        fh = open(out_vcf)
        lines = [line for line in fh.readlines() if line[0] != '#']
        linecount = len(lines)
        fh.close()

        i = 0
        # Iterate over ref names
        for ref in ('Ref1','Ref2','Ref3'):
            # Iterate over ref lengths for each refname
            for num in range(1,reflen+1):
                # Make sure all base nums are correct
                eq_([ref,str(num)], lines[i].split('\t')[:2])
                i += 1
        eq_(numrefs*reflen, linecount)

    @timed(50.0/3.0)
    @patch('ngs_mapper.base_caller.multiprocessing')
    @patch('__builtin__.open')
    @patch('ngs_mapper.base_caller.os')
    @patch('ngs_mapper.base_caller.SeqIO')
    @patch('ngs_mapper.base_caller.mpileup')
    def test_breaks_up_refs_into_chunks(self, mmpileup, mseqio, mos, mopen, mmultiprocessing):
        reflen = 100000
        numrefs = 3
        refdepth = 100
        threads = 4
        # Each reference is 1 million bases long
        ref1 = Mock(seq='A'*reflen,id='Ref1')
        ref2 = Mock(seq='T'*reflen,id='Ref2')
        ref3 = Mock(seq='G'*reflen,id='Ref3')
        mseqio.index.return_value = {'Ref1':ref1, 'Ref2': ref2, 'Ref3': ref3}
        mseqio.parse.return_value = iter([ref1, ref2, ref3])

        # Each reference is 100 reads deep
        ref1_pile = [self._mock_pileup_str('Ref1', i, 'A', refdepth, 'G'*refdepth, 'I'*refdepth, 'I'*refdepth) for i in range(1,reflen+1)]
        ref2_pile = [self._mock_pileup_str('Ref2', i, 'T', refdepth, 'T'*refdepth, 'I'*refdepth, 'I'*refdepth) for i in range(1,reflen+1)]
        ref3_pile = [self._mock_pileup_str('Ref3', i, 'G', refdepth, 'A'*refdepth, 'I'*refdepth, 'I'*refdepth) for i in range(1,reflen+1)]
        mmpileup.side_effect = [ref1_pile, ref2_pile, ref3_pile]

        with patch('ngs_mapper.base_caller.time') as time:
            out_vcf = self._C('in.bam', 'in.ref', 'out.vcf', 25, 100000, 10, 0.8, 50, 10, threads)

        expected_regionstr = [
            ('Ref1:1-25000'),
            ('Ref1:25001-50000'),
            ('Ref1:50001-75000'),
            ('Ref1:75001-100000'),
            ('Ref2:1-25000'),
            ('Ref2:25001-50000'),
            ('Ref2:50001-75000'),
            ('Ref2:75001-100000'),
            ('Ref3:1-25000'),
            ('Ref3:25001-50000'),
            ('Ref3:50001-75000'),
            ('Ref3:75001-100000'),
        ]


        i = 0
        for eregion, rcall in zip(expected_regionstr, mmultiprocessing.Process.call_args_list):
            # Call from multiprocessing.Process
            callargs = rcall[1]['args']

            # Parse out call info
            regionstr = callargs[2]
            tmpfile = callargs[3]
            
            # Every process call gets a new tempfile
            eq_('out.vcf.{0}'.format(i), tmpfile)
            eq_(eregion, regionstr)

            # Increment file count
            i += 1

class TestUnitMain(BaseInty):
    def _C( self, bamfile, reffile, vcf_output_file, regionstr=None, minbq=25, maxd=100000, mind=10, minth=0.8, biasth=50, bias=2, threads=1 ):
        from ngs_mapper.base_caller import main
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
            bias=bias,
            threads=threads
       )        
        with patch('ngs_mapper.base_caller.parse_args') as margparse:
            margparse.return_value = args
            return main()

    @attr('current')
    def test_noregion_emptypileup(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        assert_raises(InvalidRegionString, self._C, self.bam, self.ref, out_vcf, 'doesnotexist', 25, 100, 10, 0.8, 50, 2)
        #assert not self.cmp_files(self.vcf, out_vcf)

    def test_runs_multithread_noregionstring(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        r = self._C(self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8, 50, 2)
        assert self.cmp_files(self.vcf, out_vcf)

    def test_runs_single_regionstring(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')

        r = self._C(self.bam, self.ref, out_vcf+'.0', 'Ref1:2-8', 25, 100, 10, 0.8, 50, 2)
        r = self._C(self.bam, self.ref, out_vcf+'.1', 'Ref2:1-7', 25, 100, 10, 0.8, 50, 2)
        r = self._C(self.bam, self.ref, out_vcf+'.2', 'Ref3:1-8', 25, 100, 10, 0.8, 50, 2)

        # Write expected vcf by skipping some positions
        evcf = 'expected.vcf'
        with open(evcf, 'w') as fho:
            with open(self.vcf) as fhi:
                for line in fhi:
                    # Write header
                    if line.startswith('#'):
                        fho.write(line)
                        continue
                    parts = line.split('\t')
                    ref = parts[0]
                    pos = int(parts[1])

                    # Include only our test regionstr
                    if ref == 'Ref1':
                        if pos >= 2 and pos <= 8:
                            fho.write(line)
                    elif ref == 'Ref2':
                        if pos >= 1 and pos <= 7:
                            fho.write(line)
                    elif ref == 'Ref3':
                        fho.write(line)

        # Concat all generated vcf to compare with
        with open(out_vcf,'w') as fh:
            # Write all of first file
            with open(out_vcf+'.0') as fh1:
                fh.write(fh1.read())
            with open(out_vcf+'.1') as fh1:
                for line in fh1:
                    if not line.startswith('#'):
                        fh.write(line)
            with open(out_vcf+'.2') as fh1:
                for line in fh1:
                    if not line.startswith('#'):
                        fh.write(line)
        assert self.cmp_files(evcf, out_vcf)

class TestIntegrate(BaseInty):
    def _C(self, bamfile, reffile, vcf_output_file, regionstr=None, minbq=25, maxd=100000, mind=10, minth=0.8, biasth=50, bias=2, threads=2):
        import subprocess
        script_path = 'base_caller'
        cmd = [script_path, bamfile, reffile]
        if regionstr:
            cmd += ['-r', regionstr]
        if vcf_output_file:
            cmd += [vcf_output_file]
        cmd += ['-minbq', minbq, '-maxd', maxd, '-mind', mind, '-minth', minth, '-biasth', biasth, '-bias', bias, '--threads', threads]
        cmd = [str(x) for x in cmd]
        #print cmd
        return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    @timed(EXPECTED_TOTAL_TIME)
    def test_exit_0_on_success(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        p = self._C(self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8)
        eq_(0, p.wait())

    def test_outputs_correct_vcf_noregionstr(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        p = self._C(self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8)
        o,e = p.communicate()
        print o
        assert self.cmp_files(self.vcf, out_vcf)

    def test_outputs_correct_vcf_regionstr(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        p = self._C(self.bam, self.ref, out_vcf, 'Ref3:1-8', 25, 100, 10, 0.8)
        o,e = p.communicate()
        print o
        fh = open('evcf.vcf','w')
        for line in open(self.vcf):
            if line.startswith('#'):
                fh.write(line)
            else:
                _line = line.split('\t')
                if _line[0] == 'Ref3':
                    fh.write(line)
        fh.close()
        assert self.cmp_files('evcf.vcf', out_vcf)

    def test_outputs_correct_vcf_nothreads(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        p = self._C(self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8, threads=1)
        o,e = p.communicate()
        print o
        assert self.cmp_files(self.vcf, out_vcf)

    def test_stdouterr_ok_and_filesmatch(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        p = self._C(self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8)

        o,e = p.communicate()
        print "STDOUT"
        print o
        print "STDERR"
        print e
        assert e==o==''

        assert self.cmp_files(self.vcf, out_vcf)

    def test_nondefault_filesdiffer(self):
        tbam, tbai = self.temp_bam(self.bam, self.bai)
        out_vcf = join(self.tempdir, tbam + '.vcf')
        p = self._C(self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8, 50, 10)

        o,e = p.communicate()
        print "STDOUT"
        print o
        print "STDERR"
        print e
        assert e==o==''

        assert not self.cmp_files(self.vcf, out_vcf)
