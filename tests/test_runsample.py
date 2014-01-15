from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import Mock, MagicMock, patch, mock_open, call

import common
import fixtures

import os
from os.path import *
from glob import glob
import shlex
from subprocess import STDOUT, check_output, CalledProcessError

class Base(common.BaseBamRef):
    # self.bam
    # self.ref
    # self.tempdir
    reads_by_sample = ''
    f = ''
    r = ''
    @classmethod
    def setUpClass(klass):
        super(Base,klass).setUpClass()
        klass.reads_by_sample = join(klass.mytempdir,'ReadsBySample')
        os.mkdir( klass.reads_by_sample )
        f,r,ref = fixtures.get_sample_paired_reads_compressed()
        klass.f = fixtures.ungiz( f, klass.reads_by_sample )
        klass.r = fixtures.ungiz( r, klass.reads_by_sample )
        klass.ref = fixtures.ungiz( ref, klass.mytempdir )

@attr('current')
class TestUnitArgs(object):
    def _C( self, arglist ):
        from runsample import parse_args
        return parse_args( arglist )

    def test_defaults( self ):
        args = ['ReadsBySample','Reference.fasta','Sample1']
        res = self._C( args )
        eq_( 'ReadsBySample', res.readsdir )
        eq_( 'Reference.fasta', res.reference )
        eq_( 'Sample1', res.prefix )
        eq_( os.getcwd(), res.outdir )

    def test_set_outdir( self ):
        args = ['ReadsBySample','Reference.fasta','Sample1','-od','outdir']
        res = self._C( args )
        eq_( 'ReadsBySample', res.readsdir )
        eq_( 'Reference.fasta', res.reference )
        eq_( 'Sample1', res.prefix )
        eq_( 'outdir', res.outdir )

class TestFunctional(Base):
    def _run_runsample( self, readdir, reference, fileprefix, od=None ):
        script_path = dirname( dirname( abspath( __file__ ) ) )
        script_path = join( script_path, 'runsample.py' )
        cmd = script_path + ' {} {} {}'.format(readdir, reference, fileprefix)
        if od is not None:
            cmd += ' -od {}'.format(od)
        print "Running: {}".format(cmd)
        cmd = shlex.split( cmd )
        try:
            sout = check_output( cmd, stderr=STDOUT )
        except CalledProcessError as e:
            print e.output
            assert False
        return sout

    def _ensure_expected_output_files( self, outdir, prefix ):
        efiles = self._expected_files( outdir, prefix )
        for ef in efiles:
            assert isfile( ef ), "{} was not created".format(ef)
            assert os.stat( ef ).st_size > 100, "{} was less than 100 bytes".format(ef)

    def _expected_files( self, outdir, prefix ):
        #00141-98.bam  00141-98.bam.bai 00141-98.bam.qualdepth.json  00141-98.bam.qualdepth.png  00141-98.bam.qualdepth.png.pdq.tsv  00141-98.consensus.fastq  bwa.log  flagstats.txt  variants.failed.log  variants.filter.vcf  variants.indel.filter.vcf  variants.indel.raw.vcf  variants.raw.vcf
        efiles = []
        bamfile = join( outdir, prefix + '.bam' )
        efiles.append( bamfile )
        efiles.append( bamfile + '.bam.bai' )
        efiles.append( bamfile + '.bam.qualdepth.json' )
        efiles.append( bamfile + '.bam.qualdepth.png' )
        efiles.append( bamfile + '.consensus.fastq' )
        efiles.append( join( outdir, 'bwa.log' ) )
        efiles.append( join( outdir, 'flagstats.txt' ) )
        varsuffixes = ('variants.failed.log', 'variants.filter.vcf', 'variants.indel.filter.vcf', 'variants.indel.raw.vcf', 'variants.raw.vcf')
        for v in varsuffixes:
            efiles.append( join( outdir, v ) )

        return efiles

    def test_outdir_exists( self ):
        os.mkdir( 'outdir' )
        self._run_runsample( self.reads_by_sample, self.ref, 'outdir', 'tests' )
        self._ensure_expected_output_files( 'outdir', 'tests' )

    def test_outdir_not_exist( self ):
        assert not isdir( 'outdir' )
        self._run_runsample( self.reads_by_sample, self.ref, 'outdir', 'tests' )
        self._ensure_expected_output_files( 'outdir', 'tests' )
