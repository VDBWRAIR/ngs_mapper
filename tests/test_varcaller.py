from nose.tools import eq_, raises
from nose.plugins.attrib import attr
#from mock import Mock, MagicMock, patch

from os.path import *
import sys
import os
from subprocess import Popen, PIPE
import shlex
from glob import glob
import filecmp

import common
import fixtures

class Base(common.BaseClass):
    bam = join(fixtures.THIS,'fixtures','varcaller','paired.bam.gz')
    ref = join(fixtures.THIS,'fixtures','varcaller','ref.fasta.gz')
    mytempdir = ''

    @classmethod
    def setUpClass(klass):
        # Unpacks everything once so it doesn't slow down so much
        super(Base,klass).setUpClass()
        import tempfile
        klass.mytempdir = tempfile.mkdtemp()
        klass.bam = fixtures.ungiz(klass.bam,klass.mytempdir)
        klass.ref = fixtures.ungiz(klass.ref,klass.mytempdir)

    @classmethod
    def tearDownClass(klass):
        super(Base,klass).tearDownClass()
        import shutil
        shutil.rmtree(klass.mytempdir)

class TestFunctional(Base):
    expected_extensions = (
        '.failed.log', '.filter.vcf', '.indel.filter.vcf', '.indel.raw.vcf', '.raw.vcf'
    )

    def _expfiles( self, outfile ):
        return [outfile+f for f in self.expected_extensions]

    def _runsnver( self, **kwargs ):
        cmd = 'SNVer -i {infile} -r {ref} -mq {mq} -bq {bq} -s {s} -a {a} -b {b} -o {outfile}'.format(
            **kwargs
        )
        res = self._runcmd( cmd )
        print cmd
        print res[0]
        print res[1]
        return res

    def _call_varcaller( self, **kwargs ):
        d = dirname( dirname( abspath( __file__ ) ) )
        v = join( d, 'varcaller.py' )
        cmd = v + ' {infile} {ref} -mq {mq} -bq {bq} -s {s} -a {a} -b {b} -o {outfile}'.format(
            **kwargs
        )
        return self._runcmd( cmd )

    def _ensure_expected_files_exist( self, outfile ):
        # List of expected filenames
        efiles = self._expfiles( outfile )
        print "List of output files that exist"
        print glob( join( dirname( outfile ), '*' ) )
        for f in efiles:
            assert os.path.isfile( f ), "Expected output file {} was not created".format(f)

    def _cmp_expected_result_files( self, expoutfile, resoutfile ):
        self._ensure_expected_files_exist( resoutfile )
        cmp_files = zip( 
            self._expfiles( expoutfile ),
            self._expfiles( resoutfile )
        )
        for ef, rf in cmp_files:
            eq_( True, filecmp.cmp( ef, rf, 0 ), "{} and {} are not the same".format(ef,rf) )

    def _runcmd( self, cmd ):
        cmd = shlex.split( cmd )
        p = Popen( cmd, stdout=PIPE, stderr=PIPE )
        sout, serr = p.communicate()
        ret = p.returncode
        return sout, serr, ret

    def _runtst( self, infile, ref, mq, bq, s, a, b ):
        outdir = join( self.tempdir, 'varout' )
        os.mkdir( outdir )
        eoutfile = join(outdir,'expected')
        eout, eerr, eret = self._runsnver(
            infile=infile, ref=ref,
            mq=mq, bq=bq, s=s, a=a, b=b,
            outfile=eoutfile
        )

        routfile = join(outdir,'result')
        rout, rerr, rret = self._call_varcaller(
            infile=self.bam, ref=self.ref,
            mq=mq, bq=bq, s=s, a=a, b=b,
            outfile=routfile
        )
	print rout
	print rerr
        self._cmp_expected_result_files( eoutfile, routfile )
        #eq_( eout, rout )
        #eq_( eerr, rerr )
        eq_( eret, rret )

    def test_default_filters(self):
        self._runtst( self.bam, self.ref, 25, 20, 0.0001, 10, 0.2 )

    def test_modified_filters(self):
        self._runtst( self.bam, self.ref, 0, 0, 0.0001, 1, 0.05 )
