import common
from nose.tools import eq_
from mock import Mock
from os.path import *
import os
import sys
from subprocess import check_output, CalledProcessError, STDOUT
import shlex

class Base(common.BaseBamRef):
    pass

class TestUnitHandleArgs(Base):
    def _call( self, args ):
        from graphsample import handle_args
        return handle_args( args )

    def _mock_args( self ):
        args = Mock()
        args.outprefix = None
        args.outdir = self.tempdir
        args.bamfile = 'bamfile.bam'
        return args

    def test_handle_args_nothing_set(self):
        args = self._mock_args()
        res = self._call( args )
        eq_( join(self.tempdir,'bamfile.bam'), res.outpath )

    def test_handle_args_outdir_set(self):
        args = self._mock_args()
        args.outdir = 'somepath'
        res = self._call( args )
        eq_( join('somepath','bamfile.bam'), res.outpath )

    def test_handle_args_outprefix_set(self):
        args = self._mock_args()
        args.outprefix = 'someprefix'
        res = self._call( args )
        eq_( join(self.tempdir,'someprefix'), res.outpath )

    def test_handle_args_outdir_outprefix_set(self):
        args = self._mock_args()
        args.outdir = 'somepath'
        args.outprefix = 'someprefix'
        res = self._call( args )
        eq_( join('somepath','someprefix'), res.outpath )

class TestFunctional(Base):
    # Should make files with these extensions
    outfiles = ( '.qualdepth.png', '.qualdepth.json' )

    def _rungraphsample( self, bamfile, **kwargs ):
        script_path = dirname( dirname( abspath( __file__ ) ) )
        script_path = join( script_path, 'graphsample.py' )
        args = ' '.join( ['-{} {}'.format(aname,aval) for aname, aval in kwargs.items()] )
        cmd = script_path + ' {} '.format(bamfile) + args
        print "Running: {}".format(cmd)
        cmd = shlex.split( cmd )
        try:
            sout = check_output( cmd, stderr=STDOUT )
        except CalledProcessError as e:
            print e.output
            assert False
        return sout

    def _files_exist( self, outdir, outprefix ):
        out = join( outdir, outprefix )
        for of in self.outfiles:
            o = out + of
            assert os.path.isfile( o ), "Did not produce {}".format(o)

    def test_createsfiles( self ):
        res = self._rungraphsample( self.bam )
        self._files_exist( os.getcwd(), basename(self.bam) )

    def test_createsfiles_outdir_set( self ):
        os.mkdir( 'outdir' )
        res = self._rungraphsample( self.bam, od='outdir' )
        self._files_exist( 'outdir', basename(self.bam) )

    def test_createsfiles_outprefix_set( self ):
        res = self._rungraphsample( self.bam, op='outfile' )
        self._files_exist( os.getcwd(), 'outfile' )

    def test_createsfiles_outprefix_outdir_set( self ):
        os.mkdir( 'outdir' )
        res = self._rungraphsample( self.bam, od='outdir', op='outfile' )
        self._files_exist( 'outdir', 'outfile' )
