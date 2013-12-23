from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import MagicMock, patch, Mock

from common import BaseClass
import fixtures

import os
from os.path import *
from glob import glob

class Base(BaseClass):
    pass

class TestUnitBWAMem(Base):
    def test_bwa_mem_nonpaired(self):
        with patch('run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=0 ) ) ) as b:
            with patch('run_bwa.index_ref',Mock(return_value=True)) as a:
                from run_bwa import bwa_mem
                result = bwa_mem( 'F.fq', mate=None, ref='ref.fna' )
                eq_( 'bwa.sai', result )

    def test_bwa_mem_paired(self):
        with patch('run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=0 ) ) ) as b:
            with patch('run_bwa.index_ref',Mock(return_value=True)) as a:
                from run_bwa import bwa_mem
                result = bwa_mem( 'F.fq', mate='R.fq', ref='ref.fna' )
                eq_( 'bwa.sai', result )

    def test_bwa_mem_output_arg(self):
        with patch('run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=0 ) ) ) as b:
            with patch('run_bwa.index_ref',Mock(return_value=True)) as a:
                from run_bwa import bwa_mem
                result = bwa_mem( 'F.fq', mate='R.fq', ref='ref.fna', output='file.sai' )
                eq_( 'file.sai', result )

    def test_bwa_mem_fails(self):
        with patch('run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=1 ) ) ) as b:
            with patch('run_bwa.index_ref',Mock(return_value=True)) as b:
                from run_bwa import bwa_mem
                result = bwa_mem( 'F.fq', mate='R.fq', ref='ref.fna', output='file.sai' )
                eq_( 1, result )

    @patch('run_bwa.index_ref', Mock(return_value=False))
    def test_ref_index_fails(self):
        from run_bwa import bwa_mem, InvalidReference
        try:
            bwa_mem( 'F.fq', mate='R.fq', ref='ref.fna', output='file.sai' )
        except InvalidReference as e:
            pass
        else:
            assert False, "Did not raise InvalidReference"

    @patch('run_bwa.index_ref')
    @patch('run_bwa.compile_refs')
    @patch('run_bwa.BWAMem')
    def test_ref_index_directory(self, bwamem_mock, compile_refs_mock, index_ref_mock ):
        bwamem_mock.return_value.run.return_value = 1
        compile_refs_mock.return_value = 'ref_compiled.fna'
        index_ref_mock.return_value = 1
        from run_bwa import bwa_mem, InvalidReference
        ret = bwa_mem( 'F.fq', mate='R.fq', ref='ref_compiled.fna', output='file.sai' )
        eq_( 1, ret )

class TestUnitParseArgs(Base):
    def _CPA( self, argv ):
        from run_bwa import parse_args
        return parse_args( argv )

    @raises(SystemExit)
    def test_ref_reads_required( self ):
        res = self._CPA( [] )

    def test_ref_reads_set( self ):
        res = self._CPA( ['fake_read', 'fake_ref'] )
        eq_( res.reads, 'fake_read' )
        eq_( res.reference, 'fake_ref' )

    def test_platform_select_none( self ):
        res = self._CPA( ['fake_read', 'fake_ref'] )
        eq_( res.platforms, ['MiSeq','Sanger'] )

    def test_platform_select_single( self ):
        res = self._CPA( ['fake_read', 'fake_ref', '--platforms', 'Sanger'] )
        eq_( res.platforms, ['Sanger'] )

    @raises(SystemExit)
    def test_invalid_platform( self ):
        res = self._CPA( ['fake_read', 'fake_ref', '--platforms', 'invalid'] )

    def test_keep_temp_defaultoff( self ):
        res = self._CPA( ['fake_read', 'fake_ref'] )
        eq_( res.keep_temp, False )

    def test_keep_temp_set( self ):
        res = self._CPA( ['fake_read', 'fake_ref', '--keep-temp'] )
        eq_( res.keep_temp, True )

class TestIntegrateMainArgs(Base):
    pass

class TestFunctionalRunBWA(Base):
    def setUp(self):
        self.read1,self.read2,self.ref = fixtures.get_sample_paired_reads()

    def test_maps_reads_paired(self):
        from run_bwa import bwa_mem
        eq_( 'bwa.sai', bwa_mem( self.read1, self.read2, ref=self.ref ) )
        assert exists( 'bwa.sai' ), "Did not create a sai file"
        assert os.stat('bwa.sai' ).st_size != 0, "sai file created is zero bytes"

    def test_maps_reads_single(self):
        from run_bwa import bwa_mem
        eq_( 'bwa.sai', bwa_mem( self.read1, ref=self.ref ) )
        assert exists( 'bwa.sai' ), "Did not create a sai file"
        assert os.stat('bwa.sai').st_size != 0, "sai file created is zero bytes"

    def test_output_param(self):
        from run_bwa import bwa_mem
        eq_( 'file.sai', bwa_mem( self.read1, ref=self.ref, output='file.sai' ) )
        assert exists( 'file.sai' ), "Did not create a sai file"
        assert os.stat('file.sai').st_size != 0, "sai file created is zero bytes"

    def test_ref_is_directory(self):
        from run_bwa import bwa_mem
        import shutil
        os.mkdir( 'refs' )
        r1 = join('refs','ref1.fna')
        r2 = join('refs','ref2.fna')
        shutil.copy( self.ref, r1 )
        shutil.copy( self.ref, r2 )
        tot_size = os.stat(r1).st_size + os.stat(r2).st_size
        eq_( 'bwa.sai', bwa_mem( self.read1, self.read2, ref='refs' ) ) 
        # bwa.bwa.compile_refs produces reference.fa inside of current directory
        eq_( tot_size, os.stat( 'reference.fa' ).st_size )
