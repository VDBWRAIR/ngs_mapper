import random

from imports import *
from fixtures import FIXDIR

#from .. import bam_to_qualdepth as btqd
from ngs_mapper import bam_to_qualdepth as btqd

def rand_info( ):
    ''' Generates random parse_pileup() info '''
    depth = random.randint(0,1000)
    quals = [random.randint(0,40) for i in range(depth)]
    return (0,'N',depth,'N'*depth,quals)

def test_parseinfo( ):
    tst = rand_info()
    d,maxq,minq,aq = btqd.parse_info( tst )
    eq_(tst[2],d)
    eq_(max(tst[4]),maxq)
    eq_(min(tst[4]),minq)
    eq_(sum(tst[4])*1.0/len(tst[4]),aq)

def test_parseainfo( ):
    tst = [
        (0,'N',10,'N'*10,[40]*10),
        (1,'N',5,'N'*5,[20]*5),
        (2,'N',2,'N'*2,[35]*2)
    ]
    maxd,mind,maxq,minq,d,aq = btqd.parse_ainfo( tst )
    eq_(10,maxd)
    eq_(2,mind)
    eq_(40,maxq)
    eq_(20,minq)
    eq_([10,5,2],d)
    eq_([40,20,35],aq)

class Base(object):
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

class TestSetUMReads(Base):
    def setUp(self):
        super(TestSetUMReads,self).setUp()
        from ..bqd import parse_pileup
        #self.pileup = parse_pileup( mpileup( self.bamfile ) )
        self.pileup = {
            'chr1': {
                'maxd': 1000,
                'mind': 10,
                'maxq': 40,
                'minq': 10,
                'depths': [100]*10 + [1000]*10,
                'avgquals': [10]*10 + [40]*10,
                'length': 20
            }
        }

        self.idxstats = {
            'chr1': ['chr1', '20', '1000', '0'],
            'chr2': ['chr2', '20', '0', '0'],
            '*': ['*', '0', '0', '100']
        }

    def _call( self, bamfile, pileup ):
        from ..bam_to_qualdepth import set_unmapped_mapped_reads as sumr
        return sumr( bamfile, pileup )


    @patch('ngs_mapper.bam.get_refstats')
    def test_bamfile_missing_unampped_reads(self, get_refstats):
        del self.idxstats['*']
        get_refstats.return_value = self.idxstats
        res = self._call( self.bamfile, self.pileup )
        eq_(0, self.pileup['unmapped_reads'])
        eq_(1000, self.pileup['chr1']['mapped_reads'])

    @patch('ngs_mapper.bam.get_refstats')
    def test_sets_unmapped( self, get_refstats ):
        get_refstats.return_value = self.idxstats
        res = self._call( self.bamfile, self.pileup )
        eq_( 100, self.pileup['unmapped_reads'] )

    @patch('ngs_mapper.bam.get_refstats')
    def test_sets_mapped_1( self, get_refstats ):
        get_refstats.return_value = self.idxstats
        res = self._call( self.bamfile, self.pileup )
        eq_( 1000, self.pileup['chr1']['mapped_reads'] )

    @patch('ngs_mapper.bam.get_refstats')
    def test_sets_mapped_2( self, get_refstats ):
        get_refstats.return_value = self.idxstats
        res = self._call( self.bamfile, self.pileup )
        eq_( 1000, self.pileup['chr1']['mapped_reads'] )

    @patch('ngs_mapper.bam.get_refstats')
    def test_sets_reflen(self, get_refstats):
        get_refstats.return_value = self.idxstats
        res = self._call( self.bamfile, self.pileup )
        eq_( 20, self.pileup['chr1']['reflen'] )
