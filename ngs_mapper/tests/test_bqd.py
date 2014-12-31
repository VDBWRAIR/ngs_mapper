from imports import *

# Lazy import
from ngs_mapper.bqd import (
    G, N, LC, LQ, LCQ
)

class Base(BaseTester):
    modulepath = 'ngs_mapper.bqd'

    def setUp(self):
        pass

    def _make_lineargs(self, regiontypes):
        lineargs = {
            r:{'color':r,'linewidth':1} for r in regiontypes
        }
        return lineargs

    def _make_qualdepth(self, **kwargs):
        qualdepth = {
            'minq': kwargs.get('minq',0.0),
            'maxq': kwargs.get('maxq',40.0),
            'mind': kwargs.get('mind',0),
            'maxd': kwargs.get('maxd',100),
            'mapped_reads': kwargs.get('mapped_reads',1000),
            'length': kwargs.get('length',25),
            'depths': kwargs.get('depths',[
                0,0,0,0,0,      # Gap
                50,50,50,50,50, # Normal
                1,2,3,4,5,      # LowCoverage
                50,50,50,50,50, # LowQuality
                1,2,3,4,5       # LowCovQual
            ]),
            'avgquals': kwargs.get('avgquals',[
                0,0,0,0,0,      # Gap
                40,40,40,40,40, # Normal
                40,40,40,40,40, # LowCoverage
                0,5,10,15,20,   # LowQuality
                0,5,10,15,20    # LowCovQual
            ])
        }
        return qualdepth

class TestRegionsFromQualDepth(Base):
    functionname = 'regions_from_qualdepth'
        
    def setUp(self):
        super(TestRegionsFromQualDepth,self).setUp()

        # Provide a default qualdepth
        self.qualdepth = self._make_qualdepth()

    def test_returns_generator_object(self):
        import types
        r = self._C(self.qualdepth, 0, 25, 10)
        ok_( isinstance(r, types.GeneratorType), 'Did not return a Generator object' )

    def test_creates_correct_number_regions(self):
        r = list(self._C(self.qualdepth, 0, 25, 10))
        eq_(5, len(r))

    def test_creates_gap(self):
        r = list(self._C(self.qualdepth, 0, 25, 10))
        eq_(1, r[0].start)
        eq_(6, r[0].end)
        eq_(G, r[0].type)

    def test_creates_normal(self):
        r = list(self._C(self.qualdepth, 0, 25, 10))
        eq_(6, r[1].start)
        eq_(11, r[1].end)
        eq_(N, r[1].type)

    def test_creates_lowcoverage(self):
        r = list(self._C(self.qualdepth, 0, 25, 10))
        eq_(11, r[2].start)
        eq_(16, r[2].end)
        eq_(LC, r[2].type)

    def test_creates_lowquality(self):
        r = list(self._C(self.qualdepth, 0, 25, 10))
        eq_(16, r[3].start)
        eq_(21, r[3].end)
        eq_(LQ, r[3].type)

    def test_creates_lowcovqual(self):
        r = list(self._C(self.qualdepth, 0, 25, 10))
        eq_(21, r[4].start)
        eq_(26, r[4].end)
        eq_(LCQ, r[4].type)

@attr('current')
class TestGetRegionType(Base):
    functionname = 'get_region_type'

    def _check(self,t, d, q, g,dd,qq):
        eq_(t, self._C(d,q, g,dd,qq))

    def test_gap(self):
        self._check(G, 0,0, 0,25,10)
        self._check(G, 0,40, 0,25,10)

    def test_normal(self):
        self._check(N, 100,40, 0,25,10)

    def test_lowcoverage(self):
        # Check 1-9
        for i in range(1,10):
            yield self._check, LC, i,40, 0,25,10

    def test_lowqual(self):
        # Check 0-26
        for i in range(0,25):
            yield self._check, LQ, 100,i, 0,25,10

    def test_lowcovqual(self):
        for d,q in zip(range(1,10),range(0,25)):
            yield self._check, LCQ, d,q, 0,25,10

class TestLines2dFromRegion(Base):
    functionname = 'lines2d_from_regions'

    def test_returns_generator(self):
        import types
        from ngs_mapper.bqd import CoverageRegion
        ok_(isinstance(
            self._C(1, [CoverageRegion(0,0,G)], {}),
            types.GeneratorType,
        ))

    def test_creates_correct_lines(self):
        from ngs_mapper.bqd import CoverageRegion, REGIONTYPES
        regions = [
            CoverageRegion(0,0,r) for r in REGIONTYPES
        ]
        lineargs = self._make_lineargs(REGIONTYPES)
        r = self._C(0, regions, lineargs)
        for line, regiontype in zip(r,REGIONTYPES):
            eq_([0,0], line.get_ydata())
            eq_(1, line.get_linewidth())
            eq_(regiontype, line.get_color())
