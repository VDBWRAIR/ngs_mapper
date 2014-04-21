from imports import *
from datetime import date

class Base( common.BaseClass ):
    def setUp( self ):
        super( Base, self ).setUp()
        self.bc_path = join( 'Data', 'Intensities', 'BaseCalls' )

class TestGetBasecallsDir( Base ):
    def _C( self, *args, **kwargs ):
        from miseq_sync import get_basecalls_dir
        return get_basecalls_dir( *args, **kwargs )

    def test_correct_dir( self ):
        runpath = '/path/to/run/140305_M02261_0008_000000000-A6F0V'
        epath = join( runpath, self.bc_path )
        r = self._C( runpath )
        eq_( epath, r, 'Did not return correct BaseCalls path' )

class TestFileAlreadyCopied( Base ):
    def _C( self, *args, **kwargs ):
        from miseq_sync import file_already_copied
        return file_already_copied( *args, **kwargs )

    def test_same_file( self ):
        with open( 'a.txt', 'w' ) as fh:
            fh.write( 'samesize' )
        with open( 'b.txt', 'w' ) as fh:
            fh.write( 'samesize' )

        r = self._C( 'a.txt', 'b.txt' )
        ok_( r, 'Files should have been the same size' )

    def test_non_same_file( self ):
        with open( 'a.txt', 'w' ) as fh:
            fh.write( 'a.txt' )
        with open( 'b.txt', 'w' ) as fh:
            fh.write( 'b.1.txt' )

        r = self._C( 'a.txt', 'b.txt' )
        ok_( not r, 'Files were not the same, but were detected to be the same' )

class TestGetRunDate( Base ):
    def _C( self, *args, **kwargs ):
        from miseq_sync import get_rundate
        return get_rundate( *args, **kwargs )

    def test_gets_date_basename( self ):
        run = '140305_M02261_0008_000000000-A6F0V'
        r = self._C( run )
        eq_( '2014_03_05', r )

    def test_gets_date_nonbasename( self ):
        run = '/path/to/140305_M02261_0008_000000000-A6F0V'
        r = self._C( run )
        eq_( '2014_03_05', r )

class TestSamplenameFromFq( Base ):
    def _C( self, *args, **kwargs ):
        from miseq_sync import samplename_from_fq
        return samplename_from_fq( *args, **kwargs )

    def test_basenamepath( self ):
        r = self._C( 'samplename-here_S01_L001_R1_001.fastq.gz' )
        eq_( 'samplename-here', r )

    def test_non_basenamepath( self ):
        r = self._C( '/path/to/samplename-here_S01_L001_R1_001.fastq.gz' )
        eq_( 'samplename-here', r )
