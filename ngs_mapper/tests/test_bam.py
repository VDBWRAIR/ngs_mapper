from imports import *
from .fixtures import THIS, ungiz

class Base(BaseClass):
    modulepath = 'ngs_mapper.bam'

@patch('ngs_mapper.bam.subprocess.Popen')
class TestUnitMergeBam(Base):
    functionname = 'mergebams'

    samtools_cmd = ['samtools','merge','merged.bam']

    def test_input_files(self, popen_mock):
        files = ['in'+str(i) for i in range(1,3)]
        res = self._C( files, 'merged.bam' )
        eq_( [call(self.samtools_cmd+files)], popen_mock.call_args_list )

    def test_input_files_more_than_2(self, popen_mock):
        files = ['in'+str(i) for i in range(1,5)]
        res = self._C( files, 'merged.bam' )
        eq_( [call(self.samtools_cmd+files)], popen_mock.call_args_list )

    def test_input_files_lt_2(self, popen_mock):
        files = ['in1.bam']
        try:
            res = self._C( files, 'merged.bam' )
        except ValueError as e:
            pass
        else:
            assert False, "Did not raise ValueError with less than 2 bam files to merge"

    def test_input_not_list(self, popen_mock):
        files = 'in1.bam'
        try:
            res = self._C( files, 'merged.bam' )
        except ValueError as e:
            pass
        else:
            assert False, "Did not raise ValueError with non list item as argument for sortedbams"

@patch('ngs_mapper.bam.subprocess.Popen')
class TestUnitIndexBam(Base):
    functionname = 'indexbam'
    samtools_cmd = ['samtools','index']

    def test_input_existing_bam(self,popen_mock):
        res = self._C( 'sorted.bam' )
        eq_( [call(self.samtools_cmd+['sorted.bam'])], popen_mock.call_args_list )
        eq_( 'sorted.bam.bai', res )

@patch('ngs_mapper.bam.subprocess.Popen')
@patch('__builtin__.open')
class TestUnitSortBam(Base):
    functionname = 'sortbam'
    samtools_cmd = ['samtools','sort','-f','-']

    def test_input_file(self,open_mock, popen_mock):
        files = [Mock(name='file.bam'), Mock(name='file')]
        open_mock.side_effect = files
        popen_mock.return_value.stdout = None
        res = self._C( 'file.sam', 'file' )
        eq_( [call(self.samtools_cmd+['file'],stdin=files[0])], popen_mock.call_args_list )
        eq_( res, 'file' )

    def test_input_other(self, open_mock, popen_mock):
        files = [Mock(name='file')]
        open_mock.side_effect = files
        popen_mock.return_value.stdout = None
        input = Mock(spec=file)
        res = self._C( input, 'file' )
        eq_( [call(self.samtools_cmd+['file'],stdin=input)], popen_mock.call_args_list )
        eq_( res, 'file' )

    def test_output_other_fails(self, open_mock, popen_mock):
        from subprocess import PIPE
        try:
            res = self._C( 'file.sam', PIPE )
        except ValueError as e:
            pass
        else:
            assert False, "Did not raise value error on invalid output file"

@patch('ngs_mapper.bam.subprocess.Popen')
@patch('__builtin__.open')
class TestUnitSamToBam(Base):
    functionname = 'samtobam'
    samtools_cmd = ['samtools','view','-Sb','-']

    def test_input_file(self,open_mock, popen_mock):
        files = [Mock(name='file.sam'), Mock(name='file.bam')]
        open_mock.side_effect = files
        popen_mock.return_value.stdout = None
        res = self._C( 'file.sam', 'file.bam' )
        eq_( [call(self.samtools_cmd,stdin=files[0],stdout=files[1])], popen_mock.call_args_list )
        eq_( res, 'file.bam' )

    def test_input_other(self, open_mock, popen_mock):
        files = [Mock(name='file.bam')]
        open_mock.side_effect = files
        popen_mock.return_value.stdout = None
        input = Mock(spec=file)
        res = self._C( input, 'file.bam' )
        eq_( [call(self.samtools_cmd,stdin=input,stdout=files[0])], popen_mock.call_args_list )
        eq_( res, 'file.bam' )

    def test_output_other(self, open_mock, popen_mock):
        files = [Mock(name='file.sam')]
        open_mock.side_effect = files
        from subprocess import PIPE
        open_mock.return_value = 'file'
        popen_mock.return_value.stdout = PIPE
        res = self._C( 'file.sam', PIPE )
        eq_( [call(self.samtools_cmd,stdin=files[0],stdout=PIPE)], popen_mock.call_args_list )
        eq_( PIPE, res )

    def test_input_output_other(self, open_mock, popen_mock):
        from subprocess import PIPE
        open_mock.return_value = 'file'
        popen_mock.return_value.stdout = PIPE
        res = self._C( PIPE, PIPE )
        eq_( [call(self.samtools_cmd,stdin=PIPE,stdout=PIPE)], popen_mock.call_args_list )
        eq_( PIPE, res )

class TestIntegrate(Base):
    samfile = join(THIS,'fixtures','bam','samfile.sam.gz')
    unsortedbam = join(THIS,'fixtures','bam','unsorted.bam.gz')
    sortedbam = join(THIS,'fixtures','bam','sorted.bam.gz')
    bamindex = join(THIS,'fixtures','bam','sorted.bam.bai.gz')
    mergedbam = join(THIS,'fixtures','bam','merged.bam.gz')
    mergedbamindex = join(THIS,'fixtures','bam','merged.bam.bai.gz')
    mytempdir = ''

    @classmethod
    def setUpClass(klass):
        super(TestIntegrate,klass).setUpClass()
        import tempfile
        TestIntegrate.mytempdir = tempfile.mkdtemp(dir=tdir)
        TestIntegrate.samfile = ungiz(klass.samfile,klass.mytempdir)
        TestIntegrate.unsortedbam = ungiz(klass.unsortedbam,klass.mytempdir)
        TestIntegrate.sortedbam = ungiz(klass.sortedbam,klass.mytempdir)
        TestIntegrate.bamindex = ungiz(klass.bamindex,klass.mytempdir)
        TestIntegrate.mergedbam = ungiz(klass.mergedbam,klass.mytempdir)
        TestIntegrate.mergedbamindex = ungiz(klass.mergedbamindex,klass.mytempdir)

    @classmethod
    def tearDownClass(klass):
        super(TestIntegrate,klass).tearDownClass()
        import shutil
        shutil.rmtree(TestIntegrate.mytempdir)

    def _fe( self, left, right ):
        from filecmp import cmp
        assert cmp( left, right, False ), "{} was not equal to {}".format(left, right)

    def test_samtobam_pipes( self ):
        from ngs_mapper.bam import samtobam
        from subprocess import Popen, PIPE
        samfile = self.samfile
        cat = Popen(['cat',samfile],stdout=PIPE)
        pipe = samtobam( cat.stdout, PIPE )
        wc = Popen(['wc','-l'], stdin=pipe, stdout=PIPE)
        elines = len(open(self.unsortedbam).readlines())-1
        eq_( (str(elines)+'\n',None), wc.communicate() )

    def test_samtobam_files( self ):
        from ngs_mapper.bam import samtobam
        from subprocess import Popen, PIPE
        from StringIO import StringIO
        samfile = self.samfile
        bamfile = samtobam( samfile, 'output.bam' )
        eq_( 'output.bam', bamfile )
        self._fe(self.unsortedbam,'output.bam')

    def test_sortbam_pipein( self ):
        from ngs_mapper.bam import sortbam
        from subprocess import Popen, PIPE
        cat = Popen(['cat',self.unsortedbam],stdout=PIPE)
        sorted = 'expected_sorted.bam'
        bamfile = sortbam( cat.stdout, sorted )
        eq_( sorted, bamfile )
        self._fe( self.sortedbam, sorted )

    def test_sortbam_filein( self ):
        from ngs_mapper.bam import sortbam
        from subprocess import Popen, PIPE
        sorted = 'expected_sorted.bam'
        bamfile = sortbam( self.unsortedbam, sorted )
        eq_( sorted, bamfile )
        self._fe( self.sortedbam, sorted )

    def test_convert_then_sort( self ):
        from ngs_mapper.bam import samtobam, sortbam
        from subprocess import PIPE, Popen
        from gzip import open
        samfile = open(self.samfile)
        cat = Popen(['cat'], stdin=samfile, stdout=PIPE)
        convert = samtobam( cat.stdout, PIPE )
        sorted = sortbam( convert, 'sorted.bam' )
        self._fe( self.sortedbam, 'sorted.bam' )

    def test_convert_sort_index( self ):
        self.test_convert_then_sort()
        from ngs_mapper.bam import indexbam
        index = indexbam( 'sorted.bam' )
        eq_( 'sorted.bam.bai', index )
        self._fe( self.bamindex, 'sorted.bam.bai' )

    def test_convert_sort_merge_index( self ):
        self.test_convert_then_sort()
        from ngs_mapper.bam import mergebams, indexbam
        merged = mergebams( ['sorted.bam', 'sorted.bam'], 'merged.bam' )
        index = indexbam( merged )
        self._fe( self.mergedbam, 'merged.bam' )
        self._fe( self.mergedbamindex, 'merged.bam.bai' )

    def test_convert_sort_index_merge_index( self ):
        self.test_convert_sort_index()
        from ngs_mapper.bam import mergebams, indexbam
        merged = mergebams( ['sorted.bam', 'sorted.bam'], 'merged.bam' )
        index = indexbam( merged )
        self._fe( self.mergedbam, 'merged.bam' )
        self._fe( self.mergedbamindex, 'merged.bam.bai' )

### Begin better testing stuff haha
import mock
import unittest2 as unittest
import types

from .. import bam

@attr('current')
@mock.patch.object(bam, 'samtools')
@mock.patch.object(bam, 'filehandle')
@mock.patch.object(bam, 'log')
class TestBamToFastq(unittest.TestCase):
    def setUp(self):
        self.sam_rows = [
            'read1\t4\tref\t1\t60\t1171M\t=\t0\t0\tATGC\t!!!!\tAS:i:0\tXS:i:0\n',
            'read2\t4\tref\t1\t60\t1117M\t=\t0\t0\tCGTA\t!!!!\tAS:i:0\tXS:i:0\n'
        ]

    def test_missing_samtools_raises_exception(self, *args):
            mfilehandle = args[1]
            msamtools = args[2]
            mfilehandle.open.return_value = (mock.Mock(),'bam')
            excep = OSError
            msamtools.view.side_effect = excep
            self.assertRaises(
                excep,
                next, bam.bam_to_fastq('/path/to/foo.bam')
            )

    def test_accepts_filehandle(self, *args):
        msamtools = args[2]
        fh = mock.Mock(file)
        msamtools.view.side_effect = [self.sam_rows]

        r = bam.bam_to_fastq(fh)

        self.assertIsInstance(r, types.GeneratorType) 
        self.assertEqual(
            '@read1\nATGC\n+\n!!!!', next(r)
        )
        self.assertEqual(
            '@read2\nCGTA\n+\n!!!!', next(r)
        )
        msamtools.view.assert_called_once_with(fh)

    def test_accepts_filepath(self, *args):
        mfilehandle = args[1]
        msamtools = args[2]
        fh = mock.Mock(file)
        mfilehandle.open.return_value = (fh,'bam')
        f = '/path/to/foo.bam'
        msamtools.view.side_effect = [self.sam_rows]
        
        r = bam.bam_to_fastq(f)
        
        self.assertEqual(
            '@read1\nATGC\n+\n!!!!', next(r)
        )
        self.assertEqual(
            '@read2\nCGTA\n+\n!!!!', next(r)
        )

        msamtools.view.assert_called_once_with(fh)

    def test_input_contains_0_samrows_logs_warning_returns_empty_generator(self, *args):
        mfilehandle = args[1]
        msamtools = args[2]
        fh = mock.Mock(file)
        mfilehandle.open.return_value = (fh,'bam')
        mlog = args[0]
        f = '/path/to/foo.bam'
        msamtools.view.side_effect = [[]]
        
        r = bam.bam_to_fastq(f)
        self.assertIsInstance(r, types.GeneratorType) 

        self.assertRaises(
            StopIteration,
            next, r
        )

        msamtools.view.assert_called_once_with(fh)
