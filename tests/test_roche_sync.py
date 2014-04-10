from nose.tools import eq_, raises, ok_
from nose.plugins.attrib import attr

from mock import Mock, MagicMock, patch

import os
from os.path import *
import sys
import tempfile
import shutil
import tempfile
from datetime import datetime
import subprocess
import re
from glob import glob

samplesheet = '''SampleName,Region,Barcode,PrimerFileLocation
Sample1,1,Mid1,/path/to/primer.fasta
Sample2,2,RL1,'''

class BaseClass( object ):
    def setUp( self ):
        self.midparse = join(dirname(dirname(__file__)),'MidParse.conf')
        self.tdir = tempfile.mkdtemp( prefix='rochesync', suffix='tests' )
        os.chdir( self.tdir )
        self.samples = ['Sample1','Sample2']
        self.regions = ['1','2']
        self.mids = ['TI1','RL1']
        self.primers = ['/path/to/primer.fasta','']
        self.sheet = self.make_samplesheet( None, self.samples, self.regions, self.mids, self.primers )
        self.d = self.roche_dir_time()
        self.rdir = join( os.getcwd(), 'R_{}_FLX00000001_admin_00000001_test'.format(self.d) )
        self.listing = self.make_rdir( basename(self.rdir), self.samples, self.regions, self.mids, self.primers )

    def tearDown( self ):
        os.chdir( '/' )
        #shutil.rmtree( self.tdir )

    def make_filesdirs( self, root, filelist, dirlist ):
        listing = [root]
        if not isdir( root ):
            os.mkdir( root )
        for f in filelist:
            rf = join( root, f )
            open( rf, 'w' )
            listing.append( rf )
        for d in dirlist:
            rd = join( root, d )
            os.mkdir( rd )
            listing.append( rd )
        return listing

    def check_raw_synced( self, ngsdata ):
        ''' Check NGSData for synced RawData/Roche454/* '''
        raw454 = join( ngsdata, 'RawData', 'Roche454' )
        for l in self.listing:
            l = join(raw454,l)
            ok_( exists( l ), '{} did not get synced'.format(l) )

    def make_rdir( self, root, samplelist, regions, barcodes, primers ):
        '''
            Make R_ directory
            @param root - Path to R_ to create
        '''
        listing = []
        files = [
            '400x_TACG_70x75_XLPLUSKIT.icl',
            'aaLog.txt',
            'Applications_II.group',
            'calParams.tgz',
            'dataRunParams.parse',
            'imageLog.parse',
            'ptpImage.pif',
            'runComments.xml',
            'runlog.parse',
            'SampleSheet.csv',
        ]
        d = self.roche_dir_time()
        sigproc = 'D_{}_machine_signalProcessing'.format(d)
        dirs = [
            sigproc,
            'D_{}_FLX12070283_imageProcessingOnly'.format(d),
            'rawImages'
        ]
        
        listing += self.make_filesdirs( root, files, dirs )
        self.make_samplesheet( join( root, 'SampleSheet.csv' ), samplelist, regions, barcodes, primers )
        listing += self.make_sigproc( join(root, sigproc), regions )
        return listing

    def roche_dir_time( self ):
        from datetime import datetime
        return datetime.now().strftime( '%Y_%m_%d_%H_%M_%S' )

    def make_sigproc( self, root, regions ):
        '''
            Make mock signal processing dir
            @param root - Path to signal processing directory
        '''
        listing = []
        files = [
                '454AllControlMetrics.csv',
                '454BaseCallerMetrics.csv',
                '454DataProcessingDir.xml',
                '454QualityFilterMetrics.txt',
                '454RuntimeMetricsAll.txt',
                '454AllControlMetrics.txt',
                '454BaseCallerMetrics.txt',
                '454QualityFilterMetrics.csv',
                '454RuntimeMetricsAll.csv',
                'gsRunProcessor.log',
        ]
        for r in regions:
            files.append( '{}.TCA.454Reads.qual'.format(r) )
            files.append( '{}.TCA.454Reads.fna'.format(r) )
        dirs = [
            'sff',
            'regions',
        ]
        listing += self.make_filesdirs( root, files, dirs )
        listing += self.make_filesdirs( join(root,'sff'), ['ABCDEFGH0{}.sff'.format(r) for r in regions], [] )
        listing += self.make_filesdirs( join(root,'regions'), ['{}.cwf'.format(r) for r in regions], [] )
        return listing

    def make_samplesheet( self, sheetpath=None, samplelist=[], regions=[], barcodes=[], primers=[] ):
        import csv
        if sheetpath is None:
            _, sheetpath = tempfile.mkstemp( prefix='samplesheet', suffix='.csv' )
        with open( sheetpath, 'w' ) as fh:
            fh.write( 'SampleName,Region,Barcode,PrimerFileLocation\n' )
            writer = csv.writer( fh, delimiter=',' )
            for x in zip( samplelist, regions, barcodes, primers ):
                writer.writerow( x )
        return sheetpath

    def check_linkdemultiplexed( self, rdir ):
        '''
            sigprocdir needs to have R_ dir in it to pull out date
            Checks to make sure that demultiplexed directory exists inside of signalProcessing directory
            Checks to make sure that demultiplexed files got linked into ReadsBySample with correct name
            Checks to make sure ReadsBySample gets a symlink to signalProcessing
        '''
        from roche_sync import format_read_name, get_rundate, get_sigprocdir, parse_samplesheet
        ngsdatadir = dirname( dirname( dirname( rdir ) ) )
        sigprocdir = get_sigprocdir( rdir )
        demuldir = join( sigprocdir, 'demultiplexed' )
        # Ensure demultiplexed happened correctly
        self.check_demultiplexed( rdir )
        # Ensure reads by sample was populated
        rbs = join( ngsdatadir, 'ReadsBySample' )
        self.check_reads_by_sample( rdir, rbs )
        # Ensure ReadData got a symlink to signalProcessing
        ok_( isdir( join(ngsdatadir,'ReadData','Roche454',basename(sigprocdir)) ), 'Did not symlink signalProcessing into ReadData' )

    def check_demultiplexed( self, rdir ):
        '''
            Make sure sigproc has demultiplexed
            Make sure all samples have sff files
        '''
        from roche_sync import format_read_name, get_rundate, get_sigprocdir, parse_samplesheet
        d = get_rundate( rdir )
        demuldir = join( get_sigprocdir( rdir ), 'demultiplexed' )
        ok_( isdir( demuldir ), 'Missing demultiplexed directory' )
        for sample in parse_samplesheet( join(rdir,'SampleSheet.csv') ):
            sample['date'] = d
            f = format_read_name( **sample )
            f = join( demuldir, f )
            ok_( exists( f ), 'Did not create {}'.format(f) )

    def check_reads_by_sample( self, rdir, readsbysampledir ):
        from roche_sync import format_read_name, get_rundate, get_sigprocdir, parse_samplesheet
        d = get_rundate( rdir )
        sigprocdir = get_sigprocdir( rdir )
        for sample in parse_samplesheet( join(rdir,'SampleSheet.csv') ):
            samplen = sample['SampleName']
            sample['date'] = d
            f = format_read_name( **sample )
            sampledir = join( readsbysampledir, samplen )
            read = join( sampledir, f )
            readpth = join( sigprocdir, 'demultiplexed', f )
            ok_( isdir( sampledir ), 'Did not create sample directory {}'.format(sampledir) )
            # This also checks for broken symlink
            ok_( exists( read ), 'Did not create sample file {}'.format(read) )

class TestSymlinkSigProc( BaseClass ):
    def _C( self, *args, **kwargs ):
        from roche_sync import symlink_sigproc
        return symlink_sigproc( *args, **kwargs )

    def setUp( self ):
        super( TestSymlinkSigProc, self ).setUp()
        #sync( self.rdir, 
        self.ngsdata = 'NGSData'
        self.readdata = join( self.ngsdata, 'ReadData', 'Roche454' )
        os.makedirs( self.readdata )

    def test_creates_nonexisting_readdata( self ):
        from roche_sync import get_sigprocdir, sync
        os.rmdir( self.readdata )
        self._C( self.rdir, self.ngsdata )
        # Just the basename of the sigproc
        dname = basename( get_sigprocdir( self.rdir ) )
        dpath = join( 'NGSData', 'ReadData', 'Roche454', dname )
        ok_( exists( dpath ), 'Failed to create link when ReadData was missing' )

    def test_symlinks( self ):
        from roche_sync import get_sigprocdir, sync
        self._C( self.rdir, self.ngsdata )
        dname = basename( get_sigprocdir( self.rdir ) )
        lnk = join( self.readdata, dname )
        lnkdst = os.readlink( lnk )
        ok_( exists( lnk ), 'Did not link signal processing dir' )

    def test_already_exists( self ):
        self.test_symlinks()
        self.test_symlinks()

class TestSffRegionMap( BaseClass ):
    def setUp( self ):
        super( TestSffRegionMap, self ).setUp()

    def _C( self, *args, **kwargs ):
        from roche_sync import sff_region_map
        return sff_region_map( *args, **kwargs )

    def test_maps_correctly( self ):
        from roche_sync import get_sigprocdir
        sffdir = join( get_sigprocdir(self.rdir), 'sff' )
        r = self._C( sffdir )
        eq_( {'1':'ABCDEFGH01.sff', '2':'ABCDEFGH02.sff'}, r )

def d_read( insff, outsff, bname, midparse ):
    '''mock out of demultiplex_read to avoid having to have real sff files'''
    open(outsff, 'w').close()

class TestLinkReads( BaseClass ):
    def _C( self, *args, **kwargs ):
        from roche_sync import link_reads
        return link_reads( *args, **kwargs )

    def setUp( self ):
        super( TestLinkReads, self ).setUp()
        with patch('roche_sync.demultiplex_read',d_read) as dread:
            from roche_sync import demultiplex_run
            demultiplex_run( self.rdir, self.midparse )
            self.ngsdata = 'NGSData'
            self.readsbysample = join( self.ngsdata, 'ReadsBySample' )
            os.makedirs( self.readsbysample )

    def test_links_reads( self ):
        self._C( self.rdir, self.ngsdata )
        self.check_reads_by_sample( self.rdir, self.readsbysample )

    def test_already_exists( self ):
        self.test_links_reads()
        self.test_links_reads()

@patch('roche_sync.demultiplex_read',d_read)
class TestDemultiplexRun( BaseClass ):
    def setUp( self ):
        super( TestDemultiplexRun, self ).setUp()

    def _C( self, *args, **kwargs ):
        from roche_sync import demultiplex_run
        return demultiplex_run( *args, **kwargs )

    def test_demultiplexes( self ):
        from roche_sync import parse_samplesheet
        self._C( self.rdir, self.midparse )
        self.check_demultiplexed( self.rdir )

    def test_demul_dir_already_exist_not_demultiplexed( self ):
        from roche_sync import get_sigprocdir
        # Ensure if the demultiplexed directory exists, that it still creates any other files needed
        self._C( self.rdir, self.midparse )
        demuldir = join( get_sigprocdir( self.rdir ), 'demultiplexed' )
        # All sample files (filepath, stattuple)
        samplefiles = [(f,os.stat(f)) for f in glob( join(demuldir, '*') )]
        # Remove one of the samples so to test that it gets recreated
        removesample = samplefiles[0]
        noremovesample = samplefiles[1]
        os.unlink( removesample[0] )
        ok_( not exists(removesample[0]), "{} did not get removed?".format(removesample[0]) )
        # Make sure mtime will change
        import time; time.sleep(0.5);
        # Rerun demultiplex on same directory
        self._C( self.rdir, self.midparse )

        #Refetch files and stats again
        removesample_new = (removesample[0],os.stat(removesample[0]))
        noremovesample_new = (noremovesample[0],os.stat(noremovesample[0]))
        
        # Check the removed file
        eq_(
            removesample[1].st_size, removesample_new[1].st_size,
            'Size of demultiplexed removed file is not the same as it was before it was removed'
        )

        # Check the non-removed file
        eq_( 
            noremovesample[1].st_mtime, noremovesample_new[1].st_mtime,
            "Must have redemultiplexed {} as the mtime's were not the same {} -> {}".format(
                noremovesample[0], noremovesample[1].st_mtime, noremovesample_new[1].st_mtime
            )
        )

class TestParseSampleSheet( BaseClass ):
    def setUp( self ):
        super(TestParseSampleSheet, self ).setUp()

    def _C( self, *args, **kwargs ):
        from roche_sync import parse_samplesheet
        return parse_samplesheet( *args, **kwargs )

    def test_parses( self ):
        r = self._C( self.sheet )
        r1 = next(r)
        r2 = next(r)
        print r1
        eq_( self.samples[0], r1['SampleName'], 'Did not get SampleName' )
        eq_( self.regions[0], r1['Region'], 'Did not get region' )
        eq_( self.mids[0], r1['Barcode'], 'Did not get barcode' )
        eq_( self.primers[0], r1['PrimerFileLocation'], 'Did not get primer file location' )

    @raises(IOError)
    def test_missing_samplesheet( self ):
        next( self._C( '/path/not/exist.csv' ) )

class TestGetSigProcDir( BaseClass ):
    def setUp( self ):
        super(TestGetSigProcDir,self).setUp()

    def _C( self, *args, **kwargs ):
        from roche_sync import get_sigprocdir
        return get_sigprocdir( *args, **kwargs )

    def test_has_sigproc_abspath( self ):
        sig = self._C( self.rdir )
        ok_( sig.endswith('signalProcessing'), 'Returned wrong path' )
        ok_( 'R_' in dirname(sig), '{} does not contain the original rdir'.format(sig)  )

    def test_has_sigproc_basename( self ):
        sig = self._C( basename(self.rdir) )
        ok_( sig.endswith('signalProcessing'), 'Returned wrong path' )

    def test_func( self ):
        from roche_sync import sync
        sync( self.rdir, 'NGSData' )
        rdir = join( 'NGSData', 'RawData', 'Roche454', basename( self.rdir ) )
        sig = self._C( rdir )
        ok_( sig.endswith( 'signalProcessing' ), 'Returned wrong path' )

    @raises(ValueError)
    def test_error_on_missing_sigproc( self ):
        self._C( 'NGSData' )

class TestSync( BaseClass ):
    def _C( self, *args, **kwargs ):
        from roche_sync import sync
        return sync( *args, **kwargs )

    def test_syncs( self ):
        self._C( self.rdir, 'NGSData' )
        self.check_raw_synced( 'NGSData' )

    def test_skips_already_synced( self ):
        self._C( self.rdir, 'NGSData' )
        # Just check no errors
        self._C( self.rdir, 'NGSData' )

class TestFunctional( BaseClass ):
    def setUp( self ):
        super(TestFunctional,self).setUp()
        sfffile = self.mock_sfffile()
        
    def mock_sfffile( self ):
        ''' Just make fake sff files based on the -o option '''
        with open( 'sfffile', 'w' ) as fh:
            fh.write( '#!/bin/bash\n' )
            fh.write( 'touch $2\n' )
            fh.write( 'exit 0\n' )
        os.chmod( 'sfffile', 0755 )
        return 'sfffile'
        
    def _C( self, *args, **kwargs ):
        script = join( dirname( dirname(__file__) ), 'roche_sync.py' )
        cmd = 'PATH={}:$PATH {} --ngsdata {} --midparse {} {}'.format(os.getcwd(), script, args[1], args[2], args[0])
        p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True )
        eo = p.communicate()
        return (p.returncode, eo)

    def test_syncs_to_nonexisting_dst_abspath( self ):
        ngsdata = 'NGSData'
        rawdata = join( ngsdata, 'RawData', 'Roche454' )
        os.makedirs( rawdata )
        r, eo = self._C( self.rdir, ngsdata, self.midparse )
        print eo
        eq_( 0, r, 'Did not return 0' )
        self.check_raw_synced( ngsdata )
        dstrdir = join( rawdata, basename(self.rdir) )
        self.check_linkdemultiplexed( dstrdir )

    def test_syncs_to_nonexisting_basename( self ):
        ngsdata = 'NGSData'
        rawdata = join( ngsdata, 'RawData', 'Roche454' )
        os.makedirs( rawdata )
        self.rdir = basename( self.rdir )
        r, eo = self._C( self.rdir, ngsdata, self.midparse )
        print eo
        eq_( 0, r, 'Did not return 0' )
        self.check_raw_synced( ngsdata )
        dstrdir = join( rawdata, basename(self.rdir) )
        self.check_linkdemultiplexed( dstrdir )
