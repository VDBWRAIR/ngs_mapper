from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import Mock, MagicMock, patch, mock_open, call

from .common import BaseClass
from .fixtures import THIS, ungiz

import os
from os.path import *
from glob import glob

class Base(BaseClass):
    pass

@patch('bam.Popen')
@patch('__builtin__.open')
class TestFunctionalSamToBam(Base):
    samtools_cmd = ['samtools','view','-Sb','-']

    def test_input_file(self,open_mock, popen_mock):
        files = [Mock(name='file.sam'), Mock(name='file.bam')]
        open_mock.side_effect = files
        popen_mock.return_value.stdout = None
        from bam import samtobam
        res = samtobam( 'file.sam', 'file.bam' )
        eq_( [call(self.samtools_cmd,stdin=files[0],stdout=files[1])], popen_mock.call_args_list )
        eq_( res, 'file.bam' )

    def test_input_other(self, open_mock, popen_mock):
        files = [Mock(name='file.bam')]
        open_mock.side_effect = files
        popen_mock.return_value.stdout = None
        from bam import samtobam
        input = Mock(spec=file)
        res = samtobam( input, 'file.bam' )
        eq_( [call(self.samtools_cmd,stdin=input,stdout=files[0])], popen_mock.call_args_list )
        eq_( res, 'file.bam' )

    def test_output_other(self, open_mock, popen_mock):
        files = [Mock(name='file.sam')]
        open_mock.side_effect = files
        from subprocess import PIPE
        open_mock.return_value = 'file'
        popen_mock.return_value.stdout = PIPE
        from bam import samtobam
        res = samtobam( 'file.sam', PIPE )
        eq_( [call(self.samtools_cmd,stdin=files[0],stdout=PIPE)], popen_mock.call_args_list )
        eq_( PIPE, res )

    def test_input_output_other(self, open_mock, popen_mock):
        from subprocess import PIPE
        open_mock.return_value = 'file'
        popen_mock.return_value.stdout = PIPE
        from bam import samtobam
        res = samtobam( PIPE, PIPE )
        eq_( [call(self.samtools_cmd,stdin=PIPE,stdout=PIPE)], popen_mock.call_args_list )
        eq_( PIPE, res )

class TestIntegrate(Base):
    samfile = join(THIS,'fixtures','samfile.sam.gz')

    def test_samtobam_pipes( self ):
        from bam import samtobam
        from subprocess import Popen, PIPE
        samfile = ungiz(self.samfile,os.getcwd())
        cat = Popen(['cat',samfile],stdout=PIPE)
        pipe = samtobam( cat.stdout, PIPE )
        wc = Popen(['wc','-l'], stdin=pipe, stdout=PIPE)
        eq_( ('540\n',None), wc.communicate() )

    def test_samtobam_files( self ):
        from bam import samtobam
        from subprocess import Popen, PIPE
        from StringIO import StringIO
        samfile = ungiz(self.samfile,os.getcwd())
        bamfile = samtobam( samfile, 'output.bam' )
        eq_( 'output.bam', bamfile )
        assert os.path.exists('output.bam')
        assert os.stat('output.bam').st_size == 139782 
