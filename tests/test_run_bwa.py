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

class TestFunctionalRunBWA(Base):
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

class TestIntegrateRunBWA(Base):
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
