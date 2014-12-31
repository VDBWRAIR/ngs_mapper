import common
from nose.tools import eq_, raises

class TestUnitAmbTable(common.BaseClass):
    modulepath = 'ngs_mapper.alphabet'
    functionname = 'iupac_amb'

    def test_permutations( self ):
        from ngs_mapper.alphabet import AMBIGUITY_TABLE
        import itertools
        for k, v in AMBIGUITY_TABLE.items():
            for t in itertools.permutations( k ):
                yield self.tst_permutation, t, v

    def tst_permutation( self, t, v ):
        r = self._C( t )
        eq_( v, r )

    @raises(ValueError)
    def test_missing_base_raises_error( self ):
        self._C( '?' )
