from ngs_mapper.nfilter import make_filtered, write_filtered, fqs_excluding_indices, write_post_filter, mkdir_p, run_from_config

import os
import mock
import shutil
import warnings
try:
        import unittest2 as unittest
except ImportError:
        import unittest
from functools import partial
from imports import fixtures, join

class TestNGSFilter(unittest.TestCase):

    def assertListEqual(self, L1, L2):
        assert len(L1) == len(L2) and sorted(L1) == sorted(L2), "%s != %s" % (L1, L2)
#

    def setUp(self):
        fixpath = join(fixtures.THIS,'fixtures')
        fix = partial(join, fixpath)
        self.actualfn = 'testoutput/filtered.1900_S118_L001_R2_001_2015_04_24.fastq'
        self.expectedfn = fix('expected__R2__.fastq')
        self.inputfn = fix('1900/1900_S118_L001_R2_001_2015_04_24.fastq')
        self.inputdir = fix('1900')
        self.statsfile = 'testoutput/ngs_filter_stats.txt'
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
        write_post_filter(self.inputdir, 32, True, ['Sanger'], self.outdir)
        actual = open(self.actualfn)
        expected = open(self.expectedfn)
        self.assertFilesEqual(expected, actual)

    @mock.patch('ngs_mapper.nfilter.os.listdir')
    def test_fqs_excluding_indices_extensions(self, mlistdir):
        mlistdir.return_value = candidates = ['foo.fq', 'foo.sff', 'foo.fastq', 'foo.bad']
        actual = fqs_excluding_indices('D')
        expected = ['D/foo.fq', 'D/foo.fastq']
        self.assertListEqual(expected, actual)

    @mock.patch('ngs_mapper.nfilter.os.listdir')
    def test_fqs_excluding_indices_excludes_index(self, mlistdir):
        mlistdir.return_value = ['foo_1R_.fq', 'foo.sff', 'foo_I2_.fastq', 'foo__I1__.fq']
        actual = fqs_excluding_indices('D')
        expected = ['D/foo_1R_.fq' ]
        self.assertListEqual(expected, actual)

    def test_filter_raises_error_on_empty_filtered_result(self):
        ''' This should raise an AssertionError because no reads will be left after that quality filter.'''

        with warnings.catch_warnings(record=True) as w:
                # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            write_filtered(self.inputfn, 65, True)
            self.assertEquals(len(w), 1)
        #with self.assertRaises(ValueError):

    def test_symlink_file_none_filtered(self):
        write_filtered(self.inputfn, 0, False, outdir=self.outdir)
        self.assertTrue(os.path.islink(self.actualfn))

    def test_stat_file_two_filtered(self):
        write_post_filter(self.inputdir, 32, True, ['Sanger'], self.outdir)
        expected_fst = '''ngs_filter found %s reads in file %s, and filtered out %s reads.''' % (4, self.inputfn, 2)
        expected_snd = '''In file {0}, {1} reads were filtered for poor quality index below {2}. {3} reads had Ns and were filtered.'''.format(self.inputfn, 1, 32, 1)
        fst, snd =  map(str.strip, open(self.statsfile).readlines())
        self.assertEquals(fst, expected_fst)
        self.assertEquals(snd, expected_snd)

    def test_skips_platforms(self):
        try:
            write_post_filter(self.inputdir, 32, True, ['miseq'], self.outdir)
            assert False, "failed to raise valueerror on bad input"
        except ValueError, e:
            pass


    @mock.patch('ngs_mapper.nfilter.load_config')
    def test_with_config(self, mload_config):
        mfig = { 'ngs_filter' :
                {'platforms' : { 'default' : ['Sanger'] },
                 'dropNs' : { 'default' : True },
                 'indexQualityMin' : {'default' : 32},
                 'threads' : {'default' : 2}}
         }
        mload_config.return_value = mfig
        run_from_config(self.inputdir,self.outdir, '_')
        actual = open(self.actualfn)
        expected = open(self.expectedfn)
        self.assertFilesEqual(expected, actual)
