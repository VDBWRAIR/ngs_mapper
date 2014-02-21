import common
import fixtures

from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import patch, Mock, MagicMock

from StringIO import StringIO
from os.path import join, dirname, basename, exists

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
			'depth': 11,
            'A': {
                'baseq': [40],
                'mapq': [60]
            },
			'G': {
                'baseq': [1, 1, 1, 1, 40],
                'mapq': [60, 60, 60, 60, 60]
            },
			'T': {
                'baseq': [40, 1, 1, 1, 1],
                'mapq': [60, 60, 60, 60, 60]
            }
        }
        r = self._C( stats, 25, 10000, 10, 0.8 )
        eq_( ('D',3), r )

    def test_removes_N_and_updates_totaldepth( self ):
        # Bug: Returned None,0 instead of A,1
        # because when N was removed from the stats2 it was not
        # updating the total depth as well
        stats = {
            'bqualsum': 50.0,
			'mqualsum': 660.0,
			'depth': 11,
			'C': {
                'baseq': [40],
                'mapq': [60]
            },
			'T': {
                'baseq': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                'mapq': [60, 60, 60, 60, 60, 60, 60, 60, 60, 60]
            }
        }
        r = self._C( stats, 25, 100000, 10, 0.8 )
        eq_( ('C',1), r )

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

def eqs_( v1, v2, msg=None ):
    ''' Just run str on v1 and v2 and compare and use eq then '''
    if msg:
        eq_( str(v1), str(v2), msg )
    else:
        eq_( str(v1), str(v2) )

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
        eqs_( 100, r.INFO['DP'] )

    def test_depth_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eqs_( 100, r.INFO['DP'] )

    def test_refstats_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eqs_( 10, r.INFO['RC'] )
        eqs_( 10, r.INFO['PRC'] )
        eqs_( 40, r.INFO['RAQ'] )

    def test_calledbase_amb( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*9 },
            'C': { 'baseq': [38]*21 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        with patch( 'base_caller.info_stats' ) as info_stats:
            info_stats.return_value = {
                'AC':[60,21,9],
                'PAC':[60,21,10],
                'AAQ':[40,38,37],
                'bases': ['G','C','T']
            }
            r = self._C( '', 'ref:1-1', 'A', 25, 1000, 100, 0.8 )
        eq_( [60,21,9], r.INFO['AC'] )
        eq_( [60,21,10], r.INFO['PAC'] )
        eq_( [40,38,37], r.INFO['AAQ'] )
        eq_( ['G','C','T'], r.ALT )
        eq_( 'S', r.INFO['CB'] )

    def test_alternatestats_same_order( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*60 },
            'A': { 'baseq': [39]*10 },
            'C': { 'baseq': [38]*20 },
            'T': { 'baseq': [37]*10 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        with patch( 'base_caller.info_stats' ) as info_stats:
            info_stats.return_value = {
                'AC':[60,20,10],
                'PAC':[60,20,10],
                'AAQ':[40,38,37],
                'bases': ['G','C','T']
            }
            r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( [60,20,10], r.INFO['AC'] )
        eq_( [60,20,10], r.INFO['PAC'] )
        eq_( [40,38,37], r.INFO['AAQ'] )
        eq_( ['G','C','T'], r.ALT )

    def test_ensure_values_rounded( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [38]*17+[39]*16 },
            'A': { 'baseq': [30]*33 },
            'T': { 'baseq': [33]*16 + [35]*17 } # Hopefully avg qual will be 33.33 and depth is 33% too
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        r = self._C( '', 'ref:1-1', 'T', 25, 100, 10, 0.8 )
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

    def test_alternatestats_single_set( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*90 },
            'A': { 'baseq': [40]*10 },
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( [90], r.INFO['AC'] )
        eq_( [90], r.INFO['PAC'] )
        eq_( [40], r.INFO['AAQ'] )

    def test_calledbase_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 'G', r.INFO['CB'] )

    def test_calledbasedepth_set( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eqs_( 70, r.INFO['CBD'] )

    def test_fields_set_multiple( self, mock_stats ):
        mock_stats.return_value = self.mock_stats()
        # So we can check the ordering is correct
        info_stats = {'AC':[1,2,3],'AAQ':[4,5,6],'PAC':[7,8,9],'bases':['A','C','G']}
        with patch('base_caller.info_stats') as info_stats:
            info_stats.return_value = info_stats
            r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 'ref', r.CHROM )
        eqs_( 1, r.POS )
        eqs_( 'A', r.REF )
        eq_( info_stats['bases'], r.ALT )

    def test_fields_set_single( self, mock_stats ):
        base_stats = {
            'G': { 'baseq': [40]*80 },
            'A': { 'baseq': [40]*10 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        r = self._C( '', 'ref:1-1', 'A', 25, 1000, 10, 0.8 )
        eq_( 'ref', r.CHROM )
        eqs_( 1, r.POS )
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
    def setUp( self ):
        super( BaseInty, self ).setUp()
        fixpath = join( fixtures.THIS, 'fixtures', 'base_caller' )
        self.bam = join( fixpath, 'test.bam' )
        self.bai = join( fixpath, 'test.bam.bai' )
        self.ref = join( fixpath, 'testref.fasta' )
        self.sam = join( fixpath, 'test.sam' )
        self.vcf = join( fixpath, 'test.vcf' )
        self.template = join( fixpath, 'template.vcf' )

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
            print open(f2).read()
            return False

    def temp_bam( self, bam, bai ):
        import shutil
        tbam = join( self.tempdir, basename(self.bam) )
        tbai = tbam + '.bai'
        shutil.copy( self.bam, tbam )
        shutil.copy( self.bai, tbai )
        return tbam, tbai

@attr('current')
class TestUnitGenerateVCF(BaseInty):
    def _C( self, bamfile, reffile, regionstr, vcf_output_file, minbq, maxd, mind=10, minth=10, vcf_template=None ):
        from base_caller import generate_vcf, VCF_HEAD
        if vcf_template is None:
            template = VCF_HEAD.format(basename(bamfile))
        else:
            template = vcf_template
        return generate_vcf( bamfile, reffile, regionstr, vcf_output_file, 
                minbq, maxd, mind, minth, vcf_template=template ) 

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

class TestUnitMain(BaseInty):
    def _C( self, bamfile, reffile, vcf_output_file, regionstr=None, minbq=25, maxd=100000, mind=10, minth=0.8 ):
        from base_caller import main
        args = Mock(
            bamfile=bamfile,
            reffile=reffile,
            vcf_output_file=vcf_output_file,
            regionstr=regionstr,
            minbq=minbq,
            maxd=maxd,
            mind=mind,
            minth=minth
        )        
        return main( args )

    def test_runs( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        r = self._C( self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8 )
        assert self.cmp_files( self.vcf, out_vcf )

@attr('current')
class TestIntegrate(BaseInty):
    def _C( self, bamfile, reffile, vcf_output_file, regionstr=None, minbq=25, maxd=100000, mind=10, minth=0.8 ):
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
        #print cmd
        return subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE )

    def test_outputs_correct_vcf_1( self ):
        tbam, tbai = self.temp_bam( self.bam, self.bai )
        out_vcf = join( self.tempdir, tbam + '.vcf' )
        p = self._C( self.bam, self.ref, out_vcf, None, 25, 100, 10, 0.8 )
        o,e = p.communicate()
        assert self.cmp_files( self.vcf, out_vcf )

    def test_stdouterr_ok( self ):
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

