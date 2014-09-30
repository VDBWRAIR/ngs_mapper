from imports import *
from subprocess import check_output, CalledProcessError, STDOUT
import shlex

class Base(common.BaseBamRef):
    modulepath = 'miseqpipeline.graphsample'

class TestUnitHandleArgs(Base):
    functionname = 'handle_args'

    def _mock_args( self ):
        args = Mock()
        args.outprefix = None
        args.outdir = self.tempdir
        args.bamfile = 'bamfile.bam'
        return args

    def test_handle_args_nothing_set(self):
        args = self._mock_args()
        res = self._C( args )
        eq_( join(self.tempdir,'bamfile.bam'), res.outpath )

    def test_handle_args_outdir_set(self):
        args = self._mock_args()
        args.outdir = 'somepath'
        res = self._C( args )
        eq_( join('somepath','bamfile.bam'), res.outpath )

    def test_handle_args_outprefix_set(self):
        args = self._mock_args()
        args.outprefix = 'someprefix'
        res = self._C( args )
        eq_( join(self.tempdir,'someprefix'), res.outpath )

    def test_handle_args_outdir_outprefix_set(self):
        args = self._mock_args()
        args.outdir = 'somepath'
        args.outprefix = 'someprefix'
        res = self._C( args )
        eq_( join('somepath','someprefix'), res.outpath )

class TestUnitRunMontage(Base):
    functionname = 'run_montage'

    def test_expected_outputfile( self ):
        images = ['img.'+str(i)+'.png' for i in range(5)]
        outprefix = 'img.'
        expected_outfile = 'img.png'
        eo = expected_outfile
        # Touch all the files so they exist
        for f in images:
            open(f,'w').close()
        def side_effect( *args, **kwargs ):
            f = args[0][-1]
            open( f, 'w' ).close()

        with patch('subprocess.check_call') as check_call:
            check_call.side_effect = side_effect
            images = images + [expected_outfile]
            r = self._C( *images, compress='JPEG', quality=25 )
            check_call.assert_called_with( ['montage','-compress','JPEG','-quality','25']+images )
        eq_( eo, r, "Returned file path {} != {}".format(r,eo) )
        ok_( os.stat( r ), "did not create output file" )

class TestNormalizeRef(Base):
    functionname = 'normalize_ref'

    def test_replaces_punctuation( self ):
        import string
        r = self._C( string.punctuation )
        eq_( '_'*len(string.punctuation), r )

    def test_replaces_whitespace( self ):
        import string
        r = self._C( string.whitespace )
        eq_( '_'*len(string.whitespace), r )

    def test_not_replace_other( self ):
        import string
        r = self._C( string.letters + string.digits )
        ok_( '_' not in r )

class TestFunctional(Base):
    def setUp( self ):
        import fixtures
        super( TestFunctional, self ).setUp()
        self.testbam = join( fixtures.THIS, 'fixtures', 'base_caller', 'test.bam' )

    # Should make files with these extensions
    outfiles = ( '.qualdepth.png', '.qualdepth.json' )

    def _rungraphsample( self, bamfile, **kwargs ):
        script_path = 'graphsample.py'
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
        expected_files = []
        for of in self.outfiles:
            o = out + of
            expected_files.append( o )
            assert os.path.isfile( o ), "Did not produce {}".format(o)
        return expected_files

    def test_multiple_references( self ):
        res = self._rungraphsample( self.testbam )
        print res
        for f in self._files_exist( os.getcwd(), basename(self.testbam) ):
            if f.endswith( '.png' ):
                rsize = os.stat(f).st_size
                esize = 50000
                ok_( rsize > esize, "resulting size of {0} {1} was not larger than {2}".format(
                    f, rsize, esize
                ))

    def test_createsfiles( self ):
        res = self._rungraphsample( self.bam )
        self._files_exist( os.getcwd(), basename(self.bam) )

    def test_creates_imgdir( self ):
        # Puts all the images for each ref in a qualdepth dir
        from glob import glob
        res = self._rungraphsample( self.testbam )
        ok_( exists( 'qualdepth' ) )
        eq_( 'test.bam.qualdepth.png', glob( '*.png' )[0] )

    def test_outfiles_are_expected_size( self ):
        # Now just check filesizes against known
        # to make sure the graphics are correct?
        res = self._rungraphsample( self.bam )
        for f in self._files_exist( os.getcwd(), basename(self.bam) ):
            esize = 50000
            rsize = os.stat(f).st_size
            ok_( rsize > esize, "resulting size of {0} {1} was not larger than {2}".format(
                f, rsize, esize
            ))

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

    def test_no_recreate_json( self ):
        os.mkdir( 'outdir' )
        self._rungraphsample( self.bam, od='outdir', op='res' )
        es = os.stat( join('outdir','res'+self.outfiles[1]) )
        self._rungraphsample( self.bam, od='outdir', op='res', qualdepth=join('outdir','res'+self.outfiles[1]) )
        rs = os.stat( join('outdir','res'+self.outfiles[1]) )
        eq_( es, rs )
