from ngs_mapper.nfilter import make_filtered, write_filtered, fqs_excluding_indices, write_post_filter, mkdir_p
import os
import mock
import unittest
import shutil
import warnings
#
from functools import partial
from imports import fixtures, join
class TestNGSFilter(unittest.TestCase):

    def setUp(self):
        fixpath = join(fixtures.THIS,'fixtures')
        fix = partial(join, fixpath)
        self.actualfn = 'testoutput/filtered.1900_S118_L001_R2_001_2015_04_24.fastq'
        self.expectedfn = fix('expected__R2__.fastq')
        self.inputfn = fix('1900/1900_S118_L001_R2_001_2015_04_24.fastq')
        self.inputdir = fix('1900')
        self.outdir = 'testoutput'
        mkdir_p(self.outdir)

    def tearDown(self):
        if os.path.exists(self.actualfn): os.remove(self.actualfn)
        if os.path.exists(self.outdir): shutil.rmtree(self.outdir)

    def assertFilesEqual(self, f1, f2):
        if type(f1) == str: f1, f2 = open(f1), open(f2)
        l1, l2 = f1.readlines(), f2.readlines()
        #print 'f1 %s\n\nf2 %s' % (l1, l2)
        self.assertListEqual(l1, l2)

    def test_write_filtered_unfiltered(self):
        ''' input and output fastq files should be the same if filters make everything pass '''
        write_filtered(self.inputfn, 0, False, outdir=self.outdir)
        self.assertFilesEqual(self.inputfn, self.actualfn)

    def test_filter_pool(self):
        write_post_filter(self.inputdir, 32, True, self.outdir)
        actual = open(self.actualfn)
        expected = open(self.expectedfn)
        self.assertFilesEqual(expected, actual)

    @mock.patch('bioframes.nfilter.os.listdir')
    def test_fqs_excluding_indices_extensions(self, mlistdir):
        mlistdir.return_value = candidates = ['foo.fq', 'foo.sff', 'foo.fastq', 'foo.bad']
        actual = fqs_excluding_indices('D')
        expected = ['D/foo.fq', 'D/foo.sff', 'D/foo.fastq']
        self.assertListEqual(expected, actual)

    @mock.patch('bioframes.nfilter.os.listdir')
    def test_fqs_excluding_indices_excludes_index(self, mlistdir):
        mlistdir.return_value = ['foo_1R_.fq', 'foo.sff', 'foo_I2_.fastq', 'foo__I1__.fq']
        actual = fqs_excluding_indices('D')
        expected = ['D/foo_1R_.fq', 'D/foo.sff' ]
        self.assertListEqual(expected, actual)

    def test_filter_raises_error_on_empty_filtered_result(self):
        ''' This should raise an AssertionError because no reads will be left after that quality filter.'''

        with warnings.catch_warnings(record=True) as w:
                # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            write_filtered(self.inputfn, 65, True)
            self.assertEquals(len(w), 1)
        #with self.assertRaises(ValueError):

