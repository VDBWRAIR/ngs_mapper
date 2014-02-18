import common
from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import patch, Mock, MagicMock
from os.path import join, dirname, basename
import fixtures

class Base( common.BaseClass ):
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

class TestUnitLabelN(Base):
    def setUp( self ):
        super( TestUnitLabelN, self ).setUp()
        self.mock_stats()

    def _C( self, stats, minbq ):
        from base_caller import label_N
        return label_N( stats, minbq )

    def test_minq_lt( self ):
        self.stats['A']['baseq'] = [23,24,25,26]
        r = self._C( self.stats, 25 )
        eq_( [23,24], r['N']['baseq'] )

    def test_removes_empty_bases( self ):
        self.stats['A']['baseq'] = [10]*10
        self.stats['C']['baseq'] = [10]*10
        r = self._C( self.stats, 25 )
        assert 'A' not in r, 'A was not removed even though it had all < minq baseq'
        assert 'C' not in r, 'C was not removed even though it had all < minq baseq'
        eq_( [10]*20, r['N']['baseq'] )

    def test_adds_n_single_base( self ):
        self.stats['A']['baseq'] = [10,10] + [30]*6 + [10,10]
        r = self._C( self.stats, 25 )
        eq_( [10]*4, r['N']['baseq'] )

    def test_no_n( self ):
        r = self._C( self.stats, 25 )
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
        eq_( ('A',76), self._C( stats, 25, 100000, 10, 0.8 ) )

    def test_calls_multiple_ambiguous( self ):
        # 33% A, 33% G, 33% C should end up V
        base_stats = {
            'A': { 'baseq': [25]*33 },
            'G': { 'baseq': [25]*33 },
            'C': { 'baseq': [25]*33 }
        }
        stats = self.make_stats( base_stats )
        eq_( ('V',99), self._C( stats, 25, 100000, 100, 0.8 ) )

    def test_calls_low_coverage_n( self ):
        # 79% A, 21% Low Quality G should end up 21% N
        base_stats = {
            'A': { 'baseq': [25]*78 },
            'G': { 'baseq': [24]*21 }
        }
        stats = self.make_stats( base_stats )
        eq_( ('N',21), self._C( stats, 25, 100000, 100, 0.8 ) )

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

@attr('current')
@patch('base_caller.stats')
class TestUnitGenerateVcfRow(Base):
    def _C( self, bam, regionstr, refseq, minbq, maxd, mind, minth ):
        from base_caller import generate_vcf_row
        return generate_vcf_row( bam, regionstr, refseq, minbq, maxd, mind, minth )

    def mock_stats( self ):
        base_stats = {
            'G': { 'baseq': [40]*70 },
            'A': { 'baseq': [40]*10 },
            'C': { 'baseq': [40]*10 },
            'T': { 'baseq': [40]*10 }
        }
        return self.make_stats( base_stats )

    def test_regionstr_not_1( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:3-3', 'ACGT'*100, 25, 1000, 10, 0.8 )
        eq_( 70, r.info['DP'] )

    def test_depth_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 100, r.info['DP'] )

    def test_refstats_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 10, r.info['RC'] )
        eq_( 10, r.info['PRC'] )
        eq_( 40, r.info['ARQ'] )

    def test_calledbase_amb( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*9 },
            'C': { 'baseq': [38]*21 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 100, 0.8 )
        eq_( [60,21,9], r.info['AC'] )
        eq_( [60,21,10], r.info['PAC'] )
        eq_( [40,38,37], r.info['AAQ'] )
        eq_( ['G','C','T'], r.ALT )
        eq_( 'S', r.info['CB'] )

    def test_alternatestats_same_order( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*10 },
            'C': { 'baseq': [38]*20 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( [60,20,10], r.info['AC'] )
        eq_( [60,20,10], r.info['PAC'] )
        eq_( [40,38,37], r.info['AAQ'] )
        eq_( ['G','C','T'], r.ALT )

    def test_alternatestats_multiple_set( self, mock_stats ):
        stats = self.make_stats( base_stats )
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( [70,10,10], r.info['AC'] )
        eq_( [70,10,10], r.info['PAC'] )
        eq_( [40,40,40], r.info['AAQ'] )

    def test_alternatestats_single_set( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*80 },
            'A': { 'baseq': [40]*10 },
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( [80], r.info['AC'] )
        eq_( [80], r.info['PAC'] )
        eq_( [40], r.info['AAQ'] )

    def test_calledbase_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 'G', r.info['CB'] )

    def test_calledbasedepth_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 70, r.info['CBD'] )

    def test_fields_set_multiple( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 'ref', r.CHROM )
        eq_( 1, r.POS )
        eq_( 'A', r.REF )
        eq_( ['G','C','T'], r.ALT )

    def test_fields_set_single( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*80 },
            'A': { 'baseq': [40]*10 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 'ref', r.CHROM )
        eq_( 1, r.POS )
        eq_( 'A', r.REF )
        eq_( ['G'], r.ALT )

class TestUnitGenerateVCF(Base):
    # Hard to test each thing without generating sam files and vcf manually so
    # just going to let the integration tests do it...
    pass

from base_caller import InvalidRegionString
class TestUnitParseRegionString(object):
    def _C( self, regionstr ):
        from base_caller import parse_regionstring
        return parse_regionstring( regionstr )

    @raises(InvalidRegionString)
    def test_start_gt_stop( self ):
        self._C( 'ref1:2-1' )

    @raises(InvalidRegionString)
    def test_incorrect_format( self ):
        self._C( 'sometext' )
        self._C( 'sometext:1' )
        self._C( 'sometext:1-' )
        self._C( 'sometext:1- ' )
        self._C( 'sometext:a-b' )

    def test_correct_singlebase( self ):
        r = self._C( 'ref:1-1' )
        eq_( ('ref',1,1), r )

    def test_correct_multibase( self ):
        r = self._C( 'ref:1-2' )
        eq_( ('ref',1,2), r )

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

class TestIntegrate(Base):
    def setUp( self ):
        super( TestIntegrate, self ).setUp()
        fixpath = join( fixtures.THIS, 'fixtures', 'base_caller' )
        self.bam = join( fixpath, 'test.bam' )
        self.ref = join( fixpath, 'testref.fasta' )
        self.sam = join( fixpath, 'test.sam' )
        self.vcf = join( fixpath, 'test.vcf' )
        self.template = join( fixpath, 'template.vcf' )

    def _C1( self, bamfile, reffile, regionstr, vcf_output_file, 
            minbq, maxd, vcf_template, mind=10, minth=10 ):
        from base_caller import generate_vcf
        return generate_vcf( bamfile, reffile, regionstr, vcf_output_file, 
                minbq, maxd, vcf_template, mind=10, minth=10 ) 

    def _C( self, bamfile, reffile, regionstr, vcf_output_file, 
            minbq, maxd, mind=10, minth=10 ):
        import subprocess
        script_path = join( dirname( dirname( __file__ ) ) )
        script_path = join( script_path, 'base_caller.py' )
        cmd = [script_path, bamfile, reffile]
        if regionstr:
            cmd += ['-r', regionstr]
        if vcf_output_file:
            cmd += ['-o', vcf_output_file]
        cmd += ['-minbq', minbq, '-maxd', maxd, '-mind', mind, '-minth', minth]
        cmd = [str(x) for x in cmd]
        print cmd
        return subprocess.call( cmd )

    def test_correct_vcf( self ):
        import filecmp
        out_vcf = join( self.tempdir, 'out.vcf' )
        r = self._C( self.bam, self.ref, None, out_vcf, 25, 100, 10, 0.8 )
        assert filecmp.cmp( self.vcf, out_vcf, False ), "Outputted VCF was not the same as expected VCF"

    def test_default_outfile( self ):
        import filecmp
        import shutil
        tbam = join( self.tempdir, self.bam )
        shutil.copy( self.bam, tbam )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        r = self._C( tbam, self.ref, None, None, 25, 100, 10, 0.8 )
        assert filecmp.cmp( self.vcf, out_vcf, False ), "Outputted VCF was not the same as expected VCF"
