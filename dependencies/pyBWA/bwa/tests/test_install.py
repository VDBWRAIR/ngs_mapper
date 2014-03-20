from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import patch
from os.path import *
import os

import util

from .. import install

class Base( util.Base ):
    pass

class TestWhich( Base ):
    def temp_bwa( self ):
        pth = join( self.tempdir, 'bwa' )
        open( pth, 'w' ).close()
        return pth

    def test_in_path( self ):
        env = {'PATH': self.tempdir}
        pth = self.temp_bwa()
        # Executable
        os.chmod( pth, 0755 )
        with patch( 'os.environ', env ) as env:
            eq_( pth, install.which( 'bwa' ) )

    def test_not_in_path( self ):
        env = {'PATH': '/usr/bin'}
        pth = self.temp_bwa()
        with patch( 'os.environ', env ) as env:
            eq_( None, install.which( 'bwa' ) )

    def test_in_path_not_executable( self ):
        env = {'PATH': self.tempdir}
        pth = self.temp_bwa()
        # Not executable
        os.chmod( pth, 0644 )
        print os.stat( pth )
        with patch( 'os.environ', env ) as env:
            eq_( None, install.which( 'bwa' ) )
