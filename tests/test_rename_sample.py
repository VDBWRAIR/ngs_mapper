from imports import *
import common

class Base( common.BaseClass ):
    # All possible base paths that might be seen
    basepaths = (
        '/some/path/NGSData/ReadData',
        '/some/path/NGSData/RawData',
        '../../ReadData',
        '../../RawData',
    )
    # All possible variations of names
    # Paths are the ReadData paths to the reads
    platforms = {
        'Sanger': [
            ('Run_3130xl_something','sample_1_F1_2001_01_01_Den1_Den1_0001_A01.fastq'),
        ],
        'Roche454': [
            ('R_2013_10_11_16_42_14_FLX12070283_Administrator_10112013_Den1','D_2013_10_12_12_57_22_vnode_signalProcessing/demultiplexed/2/00356_06__2__RL9__2013_10_11__Den1.sff'),
            ('R_2013_10_11_16_42_14_FLX12070283_Administrator_10112013_Den1','D_2013_10_12_12_57_22_vnode_signalProcessing/demultiplexed/00356_06__2__RL9__2013_10_11__Den1.sff'),
        ],
        'MiSeq': [
            ('131217_M02261_0005_000000000-A6F07', '06564-98_S93_L001_R1_001_2013_12_17.fastq')
        ]
    }

    def print_tempdir( self ):
        print subprocess.check_output('find '+self.tempdir+' -exec ls -ld {} \;',shell=True)

    def create_ngs( self ):
        self.ngs = 'NGSData'
        self.dirs = ('ReadsBySample', 'ReadData', 'RawData')
        for d in self.dirs:
            v = join( self.ngs, d )
            setattr(self,d,v)
            os.makedirs(v)

    def make_platform_dirs( self, platform ):
        self.create_ngs()
        pdirs = []
        for d in self.dirs:
            if 'ReadsBySample' not in d:
                v = join( self.ngs, d, platform )
                os.makedirs( v )
                pdirs.append( v )
        return pdirs

class TestRunReadPath( Base ):
    def _C( self, *args, **kwargs ):
        from rename_sample import runread_path
        return runread_path( *args, **kwargs )

    def iter_plat( self, platform ):
        for bp in self.basepaths:
            for rund, read in self.platforms[platform]:
                yield join(bp,platform), rund, read

    def test_sanger( self ):
        for bp, rund, read in self.iter_plat( 'Sanger' ):
            p = join( bp, rund, read )
            r = self._C( p, 'Sanger' )
            eq_( (rund,read), r )

    def test_roche( self ):
        for bp, rund, read in self.iter_plat( 'Roche454' ):
            p = join( bp, rund, read )
            r = self._C( p, 'Roche454' )
            eq_( (rund, read), r )

class TestResolveSymlink( Base ):
    def _C( self, *args, **kwargs ):
        from rename_sample import resolve_symlink
        return resolve_symlink( *args, **kwargs )

    @raises(OSError)
    def test_notexist( self ):
        f = 'thisfile.txt'
        ok_( not exists(f) )
        r = self._C( f )

    def test_notsymlink( self ):
        # Just return normal files
        f = 'thisfile.txt'
        touch(f)
        r = self._C( f )
        eq_( f, r )

    def test_symlink_nopath( self ):
        f = 'thisfile.txt'
        l = 'thatfile.txt'
        touch(f)
        os.symlink( f, l )
        r = self._C( l )
        eq_( f, r )

    def test_symlink_relative( self ):
        f = 'thisfile.txt'
        l = join( 'dir', f )
        os.mkdir( 'dir' )
        touch(f)
        os.symlink( relpath(f,'dir'), l )
        r = self._C( l )
        eq_( f, r )

    def test_symlink_relative2( self ):
        f = join( 'dir', 'thisfile.txt' )
        l = 'thatfile.txt'
        os.mkdir( 'dir' )
        touch(f)
        os.symlink( f, l )
        r = self._C( l )
        eq_( f, r )

    def test_updown_symlink( self ):
        f = join( 'dir', 'thisfile.txt' )
        l = join( 'dir2', 'thatfile.txt' )
        os.mkdir( 'dir' )
        os.mkdir( 'dir2' )
        touch(f)
        os.symlink( relpath( f, dirname(l) ), l )
        r = self._C( l )
        eq_( f, r )

class TestRenameFile( Base ):
    def _C( self, *args, **kwargs ):
        from rename_sample import rename_file
        return rename_file( *args, **kwargs )

    def test_symlinksymlink( self ):
        f = 'somefile.txt'
        s1 = join( 'sym1', f )
        s2 = join( 'sym2', f )
        os.mkdir( 'sym1' )
        os.mkdir( 'sym2' )
        touch(f)

        # s1 link to f
        os.symlink( join('..',f), s1 )
        # s2 link to s1
        os.symlink( join('..','sym1',f), s2 )
        ok_(exists(s1))
        ok_(exists(s2))
        ok_(exists(f))
        print subprocess.check_output('find '+self.tempdir+' -exec ls -ld {} \;',shell=True)

        # Try to rename a symlink -> symlink -> file
        self._C( s2, 'file', '1' )
    
        print subprocess.check_output('find '+self.tempdir+' -exec ls -ld {} \;',shell=True)
        ok_( exists( join('sym2','some1.txt') ), 'Did not rename sym2 file' )
        ok_( exists( join('sym1','some1.txt') ), 'Did not rename sym1 file' )
        ok_( exists( 'some1.txt' ), 'Did not rename actual file' )

    def test_regularfile( self ):
        f = 'somefile.txt'
        touch( f )
        self._C( f, 'some', 'some_stuff' )
        ok_( not exists( f ), 'Old file still exists' )
        ok_( exists( 'some_stufffile.txt' ), 'File was not renamed' )

    def test_symlink( self ):
        f = join('files','somefile.txt') 
        l = 'somefile.txt'
        os.mkdir( dirname(f) )
        touch( f )
        os.symlink( f, l )
        self._C( l, 'some', 'some1' )
        print [join(root,fi) for root,dirs,files in os.walk('.') for fi in files]

        ok_( not exists(f), f + ' still exists' )
        ok_( not exists(l), l + ' still exists' )

        ok_( exists( 'some1file.txt' ), 'Link was not renamed' )
        ok_( exists(join('files','some1file.txt')), 'Actual file was not renamed' )

@attr('current')
class TestRenameMiSeq( Base ):
    def _C( self, *args, **kwargs ):
        from rename_sample import rename_miseq
        return rename_miseq( *args, **kwargs )

    def mock_miseq( self, runname, samplename ):
        read, raw = self.make_platform_dirs( 'MiSeq' )
        rbs = join( self.ReadsBySample, samplename )
        rd = join( read, runname )
        rw = join( raw, runname )
        os.mkdir(rbs)
        os.mkdir(rd)
        os.mkdir(rw)

        rawsample = '{}_S01_L001_R1_001.fastq.gz'.format(samplename)
        readsample = rawsample.replace( '.fastq.gz', '_2001_01_01.fastq' )

        bcdir = join( rw, 'Data', 'Intensities', 'BaseCalls' )
        os.makedirs( bcdir )

        rawfile = join( bcdir, rawsample )
        fqread = join( rd, readsample )
        fqrbs = join( rbs, readsample )
        touch(rawfile)
        touch(fqread)
        os.symlink( relpath(fqread,rbs), fqrbs )

        return rawfile, fqread, fqrbs

    def test_renames( self ):
        samplename = 'sample_1'
        run = 'miseqrun'
        files = self.mock_miseq( run, samplename )
        raw, read, rbs = files
        self.print_tempdir()
        self._C( rbs, samplename, samplename+'rename', 'NGSData' )
        self.print_tempdir()
        for f in files:
            ok_( not exists(f), '{} still exists'.format(f) )
            newp = f.replace(samplename,samplename+'rename')
            ok_( exists(newp), '{} did not get created' )

class TestRenameSanger( Base ):
    def _C( self, *args, **kwargs ):
        from rename_sample import rename_sanger
        return rename_sanger( *args, **kwargs )

    def mock_sanger( self, runname, samplename ):
        # Creates the full mockup of a sanger read through the data structure
        read, raw = self.make_platform_dirs( 'Sanger' )
        rbs = join( self.ReadsBySample, samplename )
        rd = join( read, runname )
        rw = join( raw, runname )
        os.mkdir(rbs)
        os.mkdir(rd)
        os.mkdir(rw)

        sample = '{}_F1_2001_01_01_Den1_Den1_0001_A01'.format(samplename)
        rawfile = join( rw, sample+'.ab1' )
        ab1read = join( rd, sample+'.ab1' )
        fqread = join( rd, sample+'.fastq' )
        ab1rbs = join( rbs, sample+'.ab1' )
        fqrbs = join( rbs, sample+'.fastq' )
        touch(rawfile)
        touch(fqread)
        os.symlink( relpath(rawfile,rd), ab1read )
        os.symlink( relpath(ab1read,rbs), ab1rbs )
        os.symlink( relpath(fqread,rbs), fqrbs )

        return rawfile, ab1read, fqread, ab1rbs, fqrbs

    def test_renames( self ):
        samplename = 'sample_1'
        run = 'Run_3130xl_something'
        files = self.mock_sanger( run, samplename )
        rawf, ab1r, fqr, ab1rbs, fqrbs = files
        self.print_tempdir()
        self._C( fqrbs, samplename, samplename+'rename', 'NGSData' )
        self._C( ab1rbs, samplename, samplename+'rename', 'NGSData' )
        print '------------------ After --------------'
        self.print_tempdir()
        for f in files:
            # Old file should not exist anymore
            ok_( not exists(f), '{} still exists'.format(f) )
            newname = f.replace(samplename,samplename+'rename')
            # SHould rename
            ok_( exists(newname), '{} does not exist'.format(newname) )
