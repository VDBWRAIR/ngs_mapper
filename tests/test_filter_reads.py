from imports import *

class FilterBase(common.BaseClass):
    def setUp( self ):
        super(FilterBase,self).setUp()
        # Contains sff, sanger and miseq paired reads
        self.read_dir = join( THIS, 'fixtures', 'reads' )
        # Only fastq files
        self.fqs = glob( join(self.read_dir,'*.fastq') )
        # Only sff files
        self.sffs = glob( join(self.read_dir,'*.sff') )
        # All reads
        self.reads = self.fqs + self.sffs

class TestTrimReadsInDir(FilterBase):
    def setUp( self ):
        super( TestTrimReadsInDir, self ).setUp()

    def _C( self, *args, **kwargs ):
        from filter_reads import trim_reads_in_dir
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

class TestTrimRead(FilterBase):
    def setUp( self ):
        super(TestTrimRead,self).setUp()

    def _C( self, *args, **kwargs ):
        from filter_reads import trim_read
        return trim_read( *args, **kwargs )

    def test_outpath_default( self ):
        # Make sure output path default option works
        for read in self.reads:
            bn = basename(read).replace('.sff','.fastq')
            self._C( read, 20 )
            try:
                os.stat(bn)
            except OSError as e:
                ok_( False, "Did not handle default out_path {}".format(e) )

    def test_trims( self ):
        # Make sure output path and returned path are ==
        # Make sure output path exists
        # Make sure output file is smaller than input file
        for read in self.reads:
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

class TestRunCutadapt(FilterBase):
    def setUp( self ):
        super(TestRunCutadapt,self).setUp()
        self.read = self.fqs[0]

    def _C( self, *args, **kwargs ):
        from filter_reads import run_cutadapt
        return run_cutadapt( *args, **kwargs )

    def test_runs_correctly( self ):
        outfq = 'output.fastq'
        r = self._C( self.read, outfq, q=20 )
        # Make sure output is correct from stderr
        ll = len(r.splitlines())
        eq_( 14, ll, 'STDERR output was not returned correctly. Got {} lines instead of 12. Output: {}'.format(ll,r) )
        # Ensure it created the correct file name
        # Stat will freak if the file does not exist
        try:
            s = os.stat( outfq )
        except IOError as e:
            ok_( False, "Did not create correct file" )
        ok_( os.stat(self.read).st_size != s.st_size, 'No trimming happened' )

class TestIntegrate(FilterBase):
    def _C( self, *args, **kwargs ):
        script = TestIntegrate.script_path('filter_reads.py')
        return TestIntegrate.run_script( '{} {} -q {} -o {}'.format(
                script, args[0], kwargs.get('q',20), kwargs.get('o','filtered_reads')
            )
        )

    def has_files( self, dir, efiles ):
        files = set( os.listdir( dir ) )
        efiles = set(efiles)
        eq_( set([]), files-efiles, "{} did not contain exactly {}".format(dir,efiles) )

    def test_runs( self ):
        r,o = self._C( self.read_dir, q=20, o='filtered_reads' )
        # Make sure exited correctly
        eq_( 0, r )
        # Make sure the file names are same as the input files
        self.has_files( 'filtered_reads', [f.replace('.sff','.fastq') for f in os.listdir(self.read_dir)] )
