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

        # Rename forward to miseq name
        f = fixtures.ungiz( f, klass.reads_by_sample )
        klass.f = join(dirname(klass.f), 'Sample1_S01_L001_R1_001_1979_01_01.fastq' )
        os.rename( f, klass.f )

        # Rename reverse to miseq name
        r = fixtures.ungiz( r, klass.reads_by_sample )
        klass.r = join(dirname(r), 'Sample1_S01_L001_R2_001_1979_01_01.fastq' )
        os.rename( r, klass.r )

        # Rename ref to correct extension
        ref = fixtures.ungiz( ref, klass.mytempdir )
        klass.ref = join( dirname(ref), 'reference.fa' )
        os.rename( ref, klass.ref )

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

class TestUnitRunCMD(object):
    import runsample
    def _C( self, cmdstr, stdin=None, stdout=None, stderr=None, script_dir=None ):
        from runsample import run_cmd
        if script_dir is None:
            return run_cmd( cmdstr, stdin, stdout, stderr )
        else:
            return run_cmd( cmdstr, stdin, stdout, stderr, script_dir )

    def test_runs_command( self ):
        import subprocess
        res = self._C( '/bin/echo test', stdout=subprocess.PIPE )
        eq_( ('test\n',None), res.communicate() )

    def test_does_pipe( self ):
        import subprocess
        p1 = self._C( '/bin/echo test', stdout=subprocess.PIPE )
        p2 = self._C( '/bin/cat -', stdin=p1.stdout, stdout=subprocess.PIPE )
        p1.stdout.close()
        eq_( ('test\n',None), p2.communicate() )

    def test_script_dir_none( self ):
        self._C( 'echo', script_dir='' )

    def test_script_dir_somepath( self ):
        self._C( 'echo', script_dir='/bin' )

    @raises(runsample.MissingCommand)
    def test_missing_cmd_exception_caught( self ):
        self._C( 'missing.sh' )

@patch('os.path.isdir')
@patch('os.path.exists')
class TestUnitTempProjdir(object):
    def _C( self, suffix, prefix ):
        from runsample import temp_projdir
        return temp_projdir( suffix, prefix )

    def test_has_tempfs( self, exists, isdir ):
        exists.return_value = True
        isdir.return_value = True
        res = self._C( 'shmtest', 'test' )
        d, bn = split( res )
        eq_( '/dev/shm', d )
        assert bn.startswith( 'shmtest' )
        assert bn.endswith( 'test' )

    def test_nothas_tempfs( self, exists, isdir ):
        exists.return_value = False
        isdir.return_value = False
        res = self._C( 'shmtest', 'test' )
        d, bn = split( res )
        eq_( '/tmp', d )
        assert bn.startswith( 'shmtest' )
        assert bn.endswith( 'test' )

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
            return (e.output,-1)
        return sout,0

    def _ensure_expected_output_files( self, outdir, prefix ):
        efiles = self._expected_files( outdir, prefix )
        ef = set( efiles )
        rf = set( [join(outdir,f) for f in os.listdir( outdir )] )
        print "Files missing from project:"
        print ef - rf
        print "Extra files in project:"
        print rf - ef
        eq_( ef, rf )
        for ef in efiles:     
            assert isfile( ef ), "{} was not created".format(ef)
            assert os.stat( ef ).st_size > 0, "{} was not > 0 bytes".format(ef)

    def _expected_files( self, outdir, prefix ):
        #00141-98.bam  00141-98.bam.bai 00141-98.bam.qualdepth.json  00141-98.bam.qualdepth.png  00141-98.bam.qualdepth.png.pdq.tsv  00141-98.consensus.fastq  bwa.log  flagstats.txt  variants.failed.log  variants.filter.vcf  variants.indel.filter.vcf  variants.indel.raw.vcf  variants.raw.vcf
        efiles = []
        bamfile = join( outdir, prefix + '.bam' )
        efiles.append( bamfile )
        efiles.append( bamfile + '.bai' )
        efiles.append( bamfile + '.qualdepth.json' )
        efiles.append( bamfile + '.qualdepth.png' )
        efiles.append( bamfile + '.consensus.fasta' )
        efiles.append( join( outdir, 'bwa.log' ) )
        efiles.append( join( outdir, 'flagstats.txt' ) )
        efiles.append( join( outdir, prefix + '.log' ) )
        efiles.append( bamfile + '.vcf' )

        return efiles

    def test_outdir_exists_nonempty_should_skip( self ):
        os.mkdir( 'outdir' )
        res,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        eq_( 0, ret )
        assert 'AlreadyExists' not in res, "Raises exception when it should not have"
        res,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        eq_( -1, ret )
        assert 'AlreadyExists' in res, "Did not raise exception"

    def test_outdir_exists_empty( self ):
        os.mkdir( 'outdir' )
        out,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        eq_( 0, ret )
        print out
        self._ensure_expected_output_files( 'outdir', 'tests' )

    def test_outdir_not_exist( self ):
        assert not isdir( 'outdir' )
        out,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        print out
        eq_( 0, ret )
        assert isdir( 'outdir' )
        self._ensure_expected_output_files( 'outdir', 'tests' )

    def test_missing_executables_exit_1( self ):
        # Change path so tempdir is first in lookup and then duplicate
        # some of the pipeline scripts and just have them return 1
        with open(join(self.tempdir,'samtools'),'w') as fh:
            fh.write( '#!/bin/bash\nexit 1\n' )
        os.chmod( join(self.tempdir,'samtools'), 0755 )
        import subprocess
        script = join( dirname( dirname( abspath( __file__ ) ) ), 'runsample.py' )
        cmd = 'export PATH={}:$PATH; {} {} {} {} -od {}'.format( self.tempdir, script, self.reads_by_sample, self.ref, 'tests', 'outdir' )
        ret = subprocess.call( cmd, shell=True )
        assert ret != 0, "Return code was 0 even though some executables returned 1"
