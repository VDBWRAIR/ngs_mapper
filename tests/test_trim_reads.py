from imports import *

class TrimBase(common.BaseClass):
    def setUp( self ):
        super(TrimBase,self).setUp()
        # Contains sanger, miseq and 454 fastq
        self.read_dir = join( THIS, 'fixtures', 'trim_reads' )
        # Only fastq files
        self.se = [
            'sample1_F1_1979_01_01_Den2_Den2_0001_A01.fastq',
            '1121__2__TI86__2012_04_04__Den2.fastq',
        ]
        self.pe = [
            ('2952_S14_L001_R1_001_2014_06_13.fastq', '2952_S14_L001_R2_001_2014_06_13.fastq'),
        ]

        self.se = [join(self.read_dir,se) for se in self.se]
        self.pe = [[join(self.read_dir,p) for p in pe] for pe in self.pe]
        # All reads
        self.reads = self.se + self.pe

class TestTrimReadsInDir(TrimBase):
    def setUp( self ):
        super( TestTrimReadsInDir, self ).setUp()

    def _C( self, *args, **kwargs ):
        from trim_reads import trim_reads_in_dir
        return trim_reads_in_dir( *args, **kwargs )

    def test_runs_correctly( self ):
        # Where to put filtered reads
        outdir = 'filtered_reads'
        # Should be all the basenames with sff replaced with fastq
        expectedfiles = [f.replace('.sff','.fastq') for f in os.listdir(self.read_dir)]
        expectedfiles = sorted(expectedfiles)
        # Run
        self._C( self.read_dir, 20, outdir )

        # Grab result files
        resultfiles = sorted(os.listdir(outdir))

        # Debugging goodies
        es,rs = set(expectedfiles), set(resultfiles)
        print "In expected not in result"
        print es - rs
        print "In result not in expected"
        print rs - es

        # Make sure lists are same
        eq_( expectedfiles, resultfiles, 'Expected files({}) was not equal to Resulting files({})'.format(expectedfiles,resultfiles) )

class TestTrimRead(TrimBase):
    def setUp( self ):
        super(TestTrimRead,self).setUp()

    def _C( self, *args, **kwargs ):
        from trim_reads import trim_read
        return trim_read( *args, **kwargs )

    def test_outpath_default( self ):
        # Make sure output path default option works
        for read in self.se:
            bn = basename(read)
            self._C( read, 20 )
            try:
                os.stat(bn)
            except OSError as e:
                ok_( False, "Did not handle default out_path {}".format(e) )

    def test_trims( self ):
        # Make sure output path and returned path are ==
        # Make sure output path exists
        # Make sure output file is smaller than input file
        for read in self.se:
            bn = basename(read)
            r = self._C( read, 40, bn )
            eq_( bn, r, 'Given outpath({}) and returned path({}) were different'.format(bn,r) )
            es = os.stat(read)
            rs = os.stat(bn)
            ok_( not samestat( es, rs ), 'Output file and inputfile are the same file' )
            ok_( 
                es.st_size > rs.st_size, 
                'Did not seem to trim the file. Output file s.st_size({}) was not smaller than input file s.st_size({})'.format(rs.st_size,es.st_size)
            )
            ok_( isdir('trim_stats'), 'Did not create trim_stats directory' )
            trimstatsfile = join( 'trim_stats', bn + '.trim_stats' )
            ok_( exists(trimstatsfile), 'Did not create trimstats file {}'.format(trimstatsfile) )

class TestRunCutadapt(TrimBase):
    def setUp( self ):
        super(TestRunCutadapt,self).setUp()
        self.read = self.se[0]

    def _C( self, *args, **kwargs ):
        from trim_reads import run_cutadapt
        return run_cutadapt( *args, **kwargs )

    def test_runs_correctly( self ):
        outfq = 'output.fastq'
        os.mkdir( 'trim_stats' )
        outstat = join( 'trim_stats', outfq + '.trim_stats' )
        r = self._C( self.read, stats=outstat, o=outfq, q=20 )
        # Make sure output is correct from stderr
        ll = len(r.splitlines())
        #eq_( 14, ll, 'STDERR output was not returned correctly. Got {} lines instead of 12. Output: {}'.format(ll,r) )
        ok_( exists(outstat), 'Did not create {} stats file'.format(outstat) )
        # Ensure it created the correct file name
        # Stat will freak if the file does not exist
        try:
            s = os.stat( outfq )
        except IOError as e:
            ok_( False, "Did not create correct file" )
        ok_( os.stat(self.read).st_size != s.st_size, 'No trimming happened' )

class TestRunTrimmomatic(TrimBase):
    def setUp( self ):
        super(TestRunTrimmomatic,self).setUp()
        os.mkdir( 'trim_stats' )
        self.outstat = join( 'trim_stats', 'output.trim_stats' )

    def _C( self, *args, **kwargs ):
        from trim_reads import run_trimmomatic
        return run_trimmomatic( *args, **kwargs )

    def test_runs_se_correctly( self ):
        for read in self.se:
            r = self._C( 'SE', read, 'output.fq', ('LEADING',20), trimlog=self.outstat )
            # Make sure output is correct from stderr
            ll = len(r.splitlines())
            #eq_( 14, ll, 'STDERR output was not returned correctly. Got {} lines instead of 12. Output: {}'.format(ll,r) )
            ok_( exists(self.outstat), 'Did not create {} stats file'.format(self.outstat) )
            # Ensure it created the correct file name
            ok_( exists('output.fq'), "Did not create correct file" )
            ok_( os.stat(read).st_size != os.stat('output.fq').st_size, 'No trimming happened' )

    def test_detects_quality_score_read( self ):
        # Make sure that it detects sanger and sets -phred33
        sanger = self.se[0]
        shutil.copy( sanger, 'different_name.fastq' )
        sanger = 'different_name.fastq'
        from data import NoPlatformFound
        try:
            out = self._C( 'SE', sanger, 'output.fastq', ('LEADING',20), trimlog=self.outstat )
            ok_( True )
        except NoPlatformFound as e:
            ok_( False, 'Raised NoPlatformFound when it should not have' )

    def test_runs_pe_correctly( self ):
        for fread, rread in self.pe:
            ofp = 'out.forward.paired.fq'
            ofu = 'out.forward.unpaired.fq'
            orp = 'out.reverse.paired.fq'
            oru = 'out.reverse.unpaired.fq'
            r = self._C( 'PE', fread, rread, ofp, ofu, orp, oru, ('LEADING',20), trimlog=self.outstat )

class TestIntegrate(TrimBase):
    def _C( self, *args, **kwargs ):
        script = TestIntegrate.script_path('trim_reads.py')
        return TestIntegrate.run_script( '{} {} -q {} -o {}'.format(
                script, args[0], kwargs.get('q',20), kwargs.get('o','trimmed_reads')
            )
        )

    def has_files( self, dir, efiles ):
        files = set( os.listdir( dir ) )
        efiles = set(efiles)
        print "Expected files: {}".format(efiles)
        print "Result files: {}".format(files)
        eq_( set([]), files-efiles, "{} did not contain exactly {}. Difference: {}".format(dir,efiles,files-efiles) )

    def test_runs( self ):
        outdir = 'trimmed_reads'
        r,o = self._C( self.read_dir, q=20, o=outdir )
        # Make sure exited correctly
        eq_( 0, r )
        print o
        # Make sure the file names are same as the input files
        self.has_files( outdir, [f.replace('.sff','.fastq') for f in os.listdir(self.read_dir)] )
        self.has_files( 'trim_stats', [f + '.trim_stats' for f in os.listdir(self.read_dir)] )
