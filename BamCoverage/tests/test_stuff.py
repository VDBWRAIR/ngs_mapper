import common
from nose.tools import eq_, ok_, raises
from mock import patch, Mock

t = Mock()
t.LCDT = 10
t.LCAQT = 10
t.GDT = 5
t.GAQT = 5

class TestRegionType(common.Base):
    def _CRT( self, d, aq ):
        from BamCoverage.bam_to_json import region_type
        return region_type( d, aq )

    @patch('BamCoverage.bam_to_json.Thresholds',t)
    def _do( self, d, q, e ):
        eq_( e, self._CRT( d, q ) )

    def test_all(self):
        tstexp = [
            (4,11,'Gap'),
            (11,4,'Gap'),
            (4,4,'Gap'),
            (5,5,'Gap'),
            (9,11,'LowCoverage'),
            (11,9,'LowQuality'),
            (10,10,'LowCoverage'),
            (11,11,'Normal'),
        ]
        for d, q, e in tstexp:
            yield self._do, d, q, e

class TestAlignmentInfo(common.Base):
    def _CAI( self, bam, rs ):
        from BamCoverage.bam import alignment_info
        return alignment_info( bam, rs )

    @patch('ngs_mapper.bam.get_refstats')
    @patch('ngs_mapper.bqd.mpileup')
    @patch('ngs_mapper.bqd.parse_pileup')
    def test_all(self,parse_pileup,mpileup,get_refstats):
        pp = {
            'ref1': dict(
                mind=1,maxd=5,minq=32.0,maxq=36.0,length=4,depths=[1,5,3,3],avgquals=[32.0,34.0,33.0,33.0]
            ),
            'ref2': dict(
                mind=1,maxd=1,minq=32.0,maxq=32.0,length=1,depths=[1],avgquals=[32.0]
            )
        }
        rs = {
            'ref1': ['ref1','5','100','10'],
            'ref2': ['ref2','10','100','10'],
            '*': ['*','0','0','0']
        }
        parse_pileup.return_value = pp
        get_refstats.return_value = rs
        g = self._CAI( 'test', '' )

        first = next(g)
        eq_( rs['ref1'], first[0] )
        eq_( pp['ref1'], first[1] )

        second = next(g)
        eq_( rs['ref2'], second[0] )
        eq_( pp['ref2'], second[1] )
