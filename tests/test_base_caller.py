import common
from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import patch, Mock, MagicMock

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

@patch('base_caller.stats')
class TestUnitCaller(Base):
    def setUp( self ):
        super( TestUnitCaller, self ).setUp()

    def _C( self, bamfile, refstr, minbq, maxd, mind, minth ):
        from base_caller import caller
        return caller( bamfile, refstr, minbq, maxd, mind, minth )

    def test_calls_multiple_ambiguous( self, mock_stats ):
        # 33% A, 33% G, 33% C should end up V
        base_stats = {
            'A': { 'baseq': [25]*33 },
            'G': { 'baseq': [25]*33 },
            'C': { 'baseq': [25]*33 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        eq_( 'V', self._C( '', '', 25, 100000, 100, 0.8 ) )

    def test_calls_low_coverage_n( self, mock_stats ):
        # 79% A, 21% Low Quality G should end up 21% N
        base_stats = {
            'A': { 'baseq': [25]*79 },
            'G': { 'baseq': [24]*21 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        eq_( 'N', self._C( '', '', 25, 100000, 100, 0.8 ) )

    def test_calls_minth( self, mock_stats ):
        # 80% A, 20% G should be called A
        base_stats = {
            'A': { 'baseq': [40]*80 },
            'G': { 'baseq': [40]*20 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        eq_( 'A', self._C( '', '', 25, 100000, 10, 0.8 ) )

    def test_calls_specific_ambiguious( self, mock_stats ):
        # 79% A, 21% G should be called R
        base_stats = {
            'A': { 'baseq': [40]*79 },
            'G': { 'baseq': [40]*21 }
        }
        stats = self.make_stats( base_stats )
        mock_stats.return_value = stats
        eq_( 'R', self._C( '', '', 25, 100000, 10, 0.8 ) )

class TestUnitCallOnPct(Base):
    pass
