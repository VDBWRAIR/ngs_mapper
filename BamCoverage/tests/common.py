from nose.tools import eq_, raises
from fixtures import FIXTURES
from os.path import *

class Base(object):
    def setUp( self ):
        self.samplename = '00103-01'
        self.bamfile = join(FIXTURES,'jsonfix','00103-01.bam')
        self.jsonfile = join(FIXTURES,'jsonfix','00103-01.json')
        self.refstats = {
            'Den4/AY618992_1/Thailand/2001/Den4_1':
                ['Den4/AY618992_1/Thailand/2001/Den4_1','10649','147751','220'],
            '*':
                ['*','0','0','23842']
        }
