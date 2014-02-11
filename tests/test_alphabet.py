import common
from nose.tools import eq_

class TestUnitAmbTable(common.BaseClass):
    def _C( self, dna ):
        from alphabet import iupac_amb
        return iupac_amb( dna )

    def test_permutations( self ):
        from alphabet import AMBIGUITY_TABLE
        import itertools
        for k, v in AMBIGUITY_TABLE.items():
            for t in itertools.permutations( k ):
                yield self.tst_permutation, t, v

    def tst_permutation( self, t, v ):
        r = self._C( t )
        eq_( r, v )
