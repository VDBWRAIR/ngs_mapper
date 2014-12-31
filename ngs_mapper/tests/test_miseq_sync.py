from imports import *
from datetime import date

class Base( common.BaseClass ):
    modulepath = 'ngs_mapper.miseq_sync'

    def setUp( self ):
        super( Base, self ).setUp()
        self.bc_path = join( 'Data', 'Intensities', 'BaseCalls' )
        self.samplesheet = join( fixtures.THIS, 'fixtures', 'SampleSheet.csv' )

class TestParseSamplesheet( Base ):
    functionname = 'parse_samplesheet'

    def test_parses_valid_sheet( self ):
        r = self._C( self.samplesheet )
        s1 = next( r )
        eq_( '001', s1['Sample_ID'] )
        eq_( 'Sample1', s1['Sample_Name'] )
        eq_( 'Plate1', s1['Sample_Plate'] )
        eq_( 'A01', s1['Sample_Well'] )
        eq_( 'N701', s1['I7_Index_ID'] )
        eq_( 'TAAGGCGA', s1['index'] )
        eq_( 'S501', s1['I5_Index_ID'] )
        eq_( 'TAGATCGC', s1['index2'] )
        # Sample id 2 & 3 have the same samplename
        s2 = next( r )
        s3 = next( r )
        eq_( s2['Sample_Name'], s3['Sample_Name'] )
        s4 = next( r )
        s5 = next( r )
        try:
            next( r )
            ok_( False, 'Should not be 6 samples, but there were' )
        except:
            pass

class TestGetBasecallsDir( Base ):
    functionname = 'get_basecalls_dir'

    def test_correct_dir( self ):
        runpath = '/path/to/run/140305_M02261_0008_000000000-A6F0V'
        epath = join( runpath, self.bc_path )
        r = self._C( runpath )
        eq_( epath, r, 'Did not return correct BaseCalls path' )

class TestFileAlreadyCopied( Base ):
    functionname = 'file_already_copied'

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
    functionname = 'get_rundate'

    def test_gets_date_basename( self ):
        run = '140305_M02261_0008_000000000-A6F0V'
        r = self._C( run )
        eq_( '2014_03_05', r )

    def test_gets_date_nonbasename( self ):
        run = '/path/to/140305_M02261_0008_000000000-A6F0V'
        r = self._C( run )
        eq_( '2014_03_05', r )

    @raises(ValueError)
    def test_bad_has_no_rundate( self ):
        run = 'ABCDEFG'
        r = self._C( run )

class TestSamplenameFromFq( Base ):
    functionname = 'samplename_from_fq'

    def test_basenamepath( self ):
        r = self._C( 'samplename-here_S01_L001_R1_001.fastq.gz' )
        eq_( 'samplename-here', r )

    def test_non_basenamepath( self ):
        r = self._C( '/path/to/samplename-here_S01_L001_R1_001.fastq.gz' )
        eq_( 'samplename-here', r )

class FunctionalBase( Base ):
    def setUp( self ):
        super( FunctionalBase, self ).setUp()

    def mock_samples( self, samples, outdir ):
        ''' Make mock sample files in outdir from '''
        for sampleid, sampleinfo in self.samples.items():
            readpaths = self.mock_miseq_reads_for_sample( sampleinfo['Sample_Name'], sampleid, outdir )

    def mock_miseq_reads_for_sample( self, samplename, sampleid, outpath ):
        import gzip
        template = '{samplename}_S{sampleid}_L001_R{rf}_001.fastq.gz'
        reads = []
        for i in (1,2):
            readpath = join( outpath, template.format( samplename=samplename, sampleid=sampleid, rf=i ) )
            with gzip.open( readpath, 'wb' ) as fh:
                fh.write( '@{}\nATGC\n+\nIIII\n'.format(samplename) )
            reads.append( readpath )
        return reads
        
    def mock_miseq_run( self ):
        from ngs_mapper.miseq_sync import get_basecalls_dir, parse_samplesheet
        # Fixture samplesheet file
        self.rundate = '010101'
        self.runname = '{}_M02261_0001_00000000-A6F0V'.format(self.rundate)
        self.bcdir = get_basecalls_dir( self.runname )
        # Create all directories to BaseCalls dir
        os.makedirs( self.bcdir )
        # Copy in the samplesheet
        shutil.copy( self.samplesheet, self.runname )
        self.samplesheet = join( self.runname, 'SampleSheet.csv' )
        # Parse the samplesheet
        self.samples = {sample['Sample_ID']:sample for sample in parse_samplesheet(self.samplesheet)}
        self.mock_samples( self.samples, self.bcdir )

class TestFunctional( FunctionalBase ):
    def _C( self, *args, **kwargs ):
        script = 'miseq_sync.py'
        cmd = [script, '--ngsdata', args[0], args[1]]
        return subprocess.call( cmd )

    def ensure_fastqfiles( self, runpath ):
        from ngs_mapper.miseq_sync import get_basecalls_dir, parse_samplesheet
        bcdir = get_basecalls_dir( runpath )
        samples = list( parse_samplesheet( join(runpath, 'SampleSheet.csv') ) )
        for sample in samples:
            sample_fqs = glob( join(bcdir, '{}_S{}*'.format(sample['Sample_Name'],sample['Sample_ID'])) )
            eq_( 2, len(sample_fqs) )

    def ensure_samplesheet( self, runpath ):
        ok_( exists( join(runpath,'SampleSheet.csv') ), '{} missing samplesheet'.format(runpath) )

    def ensure_data_structure( self, runpath, ngsdata ):
        from ngs_mapper.miseq_sync import parse_samplesheet
        listing = [join(root,file) for root, files, dirs in os.walk(ngsdata) for file in files]
        print listing
        print open( 'pipeline.log' ).read()
        rawdata = join( ngsdata, 'RawData', 'MiSeq', basename(runpath) )
        readdata = join( ngsdata, 'ReadData', 'MiSeq', basename(runpath) )
        readsbysample = join( ngsdata, 'ReadsBySample' )
        ok_( exists(rawdata), 'Did not create {}'.format(rawdata) )
        ok_( exists(readdata), 'Did not create {}'.format(readdata) )
        # Verifies SampleSheet.csv was copied
        samples = list( parse_samplesheet( join(rawdata,'SampleSheet.csv') ) )
        for sample in samples:
            self.ensure_reads_by_sample( readsbysample, readdata, sample['Sample_Name'], sample['Sample_ID'] )

    def ensure_reads_by_sample( self, readsbysampledir, readdata, samplename, sampleid ):
        # Should be links
        rbs = join( readsbysampledir, samplename )
        ok_( exists(rbs), 'Did not create {}'.format(rbs) )
        sample_fqs = glob( join( rbs, '{}_S{}*.fastq'.format(samplename,sampleid) ) )
        eq_( 2, len(sample_fqs), 'There were not 2 reads for {}. Reads found in {}: {}'.format(samplename,rbs,sample_fqs) )
        # Every fq in readsbysample should be a link to the runname in readdata
        for fq in sample_fqs:
            lnkdst = os.readlink( fq )
            runname = basename( dirname( lnkdst ) )
            runpth = join( '..', '..', 'ReadData', 'MiSeq', runname )
            ok_( exists( fq ), '{} is a broken link to {}'.format(fq, lnkdst) )
            eq_( runpth, dirname(lnkdst), '{} is not a correct lnk path. It should be {}'.format(lnkdst, runpth) )

    def test_ensure_synced( self ):
        self.mock_miseq_run()
        ngsdata = 'NGSData'
        eq_( 0, self._C( ngsdata, self.runname ), 'miseq_sync.py did not return 0' )
        self.ensure_data_structure( self.runname, ngsdata )

    def test_fixes_issue_1060_runpath_has_slash_at_end( self ):
        self.mock_miseq_run()
        ngsdata = 'NGSData'
        runpathwithslash = join( os.getcwd(), self.runname ) + '/'
        eq_( 0, self._C( ngsdata, runpathwithslash ), 'miseq_sync.py did not return 0' )
        self.ensure_data_structure( self.runname, ngsdata )
