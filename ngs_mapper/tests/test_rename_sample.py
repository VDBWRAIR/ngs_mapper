from imports import *
import common
from ngs_mapper.rename_sample import RenameException

class Base( common.BaseClass ):
    modulepath = 'ngs_mapper.rename_sample'

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
        ],
        'IonTorrent': [
            ('2001_01_01', 'SFFCreator_out/IonXpress_017_R_2013_01_19_15_14_09_user_VDB-47_Auto_user_VDB-47_45.sff'),
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
            try:
                os.makedirs(v)
            except OSError as e:
                if not e.errno == 17:
                    raise e

    def make_platform_dirs( self, platform ):
        self.create_ngs()
        pdirs = []
        for d in self.dirs:
            if 'ReadsBySample' not in d:
                v = join( self.ngs, d, platform )
                os.makedirs( v )
                pdirs.append( v )
        return pdirs

    def make_dirs( self, platform, runname, samplename ):
        read, raw = self.make_platform_dirs( platform )
        rbs = join( self.ReadsBySample, samplename )
        rd = join( read, runname )
        rw = join( raw, runname )
        if not isdir(rbs):
            os.mkdir(rbs)
        if not isdir(rd):
            os.mkdir(rd)
        if not isdir(rw):
            os.mkdir(rw)
        return rw, rd, rbs

    def mock_sanger( self, runname, samplename ):
        rw, rd, rbs = self.make_dirs('Sanger', runname, samplename)
        
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

    def mock_iontorrent( self, runname, samplename ):
        rw, rd, rbs = self.make_dirs('IonTorrent', runname, samplename)
        rawfile = join( rw, 'somedirs', 'IonXpress_001_R_2001_01_01_01_01_01_user_XXX-01_Auto_user_XXX-01_02.sff' )
        rawsample = '{}__1__RL1__2001_01_01__Den1.sff'.format(samplename)
        samplefile = samplename + '__1__IX001__2001_01_01__vir.sff' 
        rdfile = join( rd, samplefile )
        rbsfile = join( rbs, samplefile )

        os.makedirs( dirname(rawfile) )
        touch( rawfile )

        rawfile_rel = relpath( rawfile, dirname(rdfile) )
        rdfile_rel = relpath( rdfile, dirname(rbsfile) )
        os.symlink( rawfile_rel, rdfile )
        os.symlink( rdfile_rel, rbsfile )

        return rawfile, rdfile, rbsfile

    def mock_roche( self, runname, samplename, byregion=True, rbsraw=False ):
        '''
            if byregion is True then demultiplex will have a region directory where the read is placed
            otherwise the read is placed directly inside of the demultiplexed directory

            if rbsraw is True then symlink rbsfile directly to RawData otherwise symlink to readdata
        '''
        rw, rd, rbs = self.make_dirs('Roche454', runname, samplename)
        # 454 is different since the R_ directory does not exist in ReadData, but instead
        # the *signalProcessing does as a symlink to RawData
        os.rmdir( rd )
        ddir = 'D_2001_01_01_01_01_01_vnode_signalProcessing'
        # rd is just the sigproc dir not the R_ dir so we fix it here
        rd = join( self.ReadData, 'Roche454', ddir )

        rawsample = '{}__1__RL1__2001_01_01__Den1.sff'.format(samplename)

        if byregion:
            rwdemul = join( rw, ddir, 'demultiplexed', '1' )
            rddemul = join( rd, 'demultiplexed', '1' )
        else:
            rwdemul = join( rw, ddir, 'demultiplexed' )
            rddemul = join( rd, 'demultiplexed' )
        # Create demultiplexed directory in RawData
        os.makedirs( rwdemul )

        # Do the sigproc symlink
        target = relpath( join(rw,ddir), join(self.ReadData,'Roche454') )
        os.symlink( target, rd )
        ok_( exists(rd), 'Sigproc {} -> {} symlink broken'.format(rd, target) )

        ok_( exists(rwdemul), '{} does not exist'.format(rwdemul) )
        ok_( exists(rddemul), '{} does not exist'.format(rddemul) )

        rawfile = join( rwdemul, rawsample )
        touch(rawfile)

        readfile = join( rddemul, rawsample )
        ok_( exists(dirname(readfile)),  )
        ok_( exists(readfile) )
        ok_( samefile(rawfile,readfile) )

        rbsfile = join( rbs, rawsample )
        if rbsraw:
            os.symlink( relpath(rawfile,rbs), rbsfile )
        else:
            os.symlink( relpath(readfile,rbs), rbsfile )
        ok_( exists(rbsfile) )

        return rawfile, readfile, rbsfile

    def mock_miseq( self, runname, samplename ):
        rw, rd, rbs = self.make_dirs('MiSeq', runname, samplename)

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


class TestRunReadPath( Base ):
    functionname = 'runread_path'

    def iter_plat( self, platform ):
        for bp in self.basepaths:
            for rund, read in self.platforms[platform]:
                yield join(bp,platform), rund, read

    def test_iontorrent( self ):
        for bp, rund, read in self.iter_plat( 'IonTorrent' ):
            p = join( bp, rund, read )
            r = self._C( p, 'IonTorrent' )
            eq_( (rund,read), r )

    def test_miseq( self ):
        for bp, rund, read in self.iter_plat( 'MiSeq' ):
            p = join( bp, rund, read )
            r = self._C( p, 'MiSeq' )
            eq_( (rund,read), r )

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
    functionname = 'resolve_symlink'

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
    functionname = 'rename_file'

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
        self._C( s2, 'somefile', 'some1' )
    
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

    def test_dir_contains_fromstr( self ):
        # Only replace when the replace string is immediately preceeded by os.sep
        f = join( 'Run_3130xl', '313', '313.txt' )
        os.makedirs( dirname(f) )
        touch( f )
        e = join( 'Run_3130xl', 'replaced', 'replaced.txt' )
        os.makedirs(dirname(e))
        self._C( f, '313', 'replaced' )
        ok_( not exists(f.replace('313','replaced')), 'Incorrectly replaced replace string' )
        ok_( exists(e), '{} was not created'.format(e) )

    def test_missing_originalname( self ):
        # Sometimes the file that the symlinks resolve to doesn't contain the
        # fromname. That file can be left alone
        os.mkdir( 'RBS' )
        os.mkdir( 'ReadD' )
        os.mkdir( 'RawD' )
        # Make the dirs
        rawf = join( 'RawD', 'filename.sff' )
        readf = join( 'ReadD', 'fromname.sff' )
        rbsf = join( 'RBS', 'fromname.sff' )
        # Make the files
        touch( rawf )
        os.symlink( relpath(rawf, 'ReadD'), readf )
        os.symlink( relpath(readf, 'RBS'), rbsf )

        self._C( rbsf, 'fromname', 'tothis' )

        self.print_tempdir()
        ok_( not exists(rbsf), 'Did not rename {}'.format(rbsf) )
        ok_( exists(rbsf.replace('fromname','tothis')), 'Did not create new symlink' )
        ok_( not exists(readf), 'Did not rename intermediate symlink file' )
        ok_( exists(readf.replace('fromname','tothis')), 'Did not recreate intermediate symlink file' )
        ok_( exists(rawf), 'Actual file no longer exiss' )

    @raises(RenameException)
    def test_newp_exists( self ):
        # Existing file
        f = join( 'Run_3130xl', 'file1.txt' )
        ref = join( 'file1', 'file1.txt' )
        os.makedirs( dirname(f) )
        os.makedirs( dirname(ref) )
        touch( f )
        os.symlink( f, ref )
        # File to rename to existing file
        f2 = join( 'Run_3130xl', 'file2.txt' )
        f2s = join( 'file2', 'file2.txt' )
        os.makedirs( dirname(f2s) )
        touch( f2 )
        os.symlink( f2, f2s )
        self._C( f2, 'file2', 'file1' )
        self.print_tempdir()
        ok_( exists(f2s), 'Broken symlink' )
        ok_( exists(f2), 'Removed f2 when it should not have' )
        ok_( samefile(f,f2), 'Overwrote f with f2' )

class TestRenameIonTorrent( Base ):
    functionname = 'rename_iontorrent'

    def test_renames_correctly( self ):
        raw, read, rbs = self.mock_iontorrent( 'somerunname', 'fromsamplename' ) 

        self._C( rbs, 'fromsamplename', 'tosamplename', 'NGSData' )

        newraw = raw
        newread = read.replace('fromsamplename', 'tosamplename')
        newrbs = rbs.replace('fromsamplename', 'tosamplename')

        self.print_tempdir()
        ok_( not exists(read), 'Did not rename the ReadData symlink({} still exists)'.format(read) )
        ok_( not exists(rbs), 'Did not rename ReadsBySample symlink({} still exists)'.format(rbs) )

        ok_( exists(raw), 'Renamed IonXpress raw file when it should not have' )
        ok_( exists(newread), 'Did not create the correct ReadData symlink' )
        ok_( exists(newrbs), 'Did not create the correct ReadsBySample symlink' )

class TestRenameRoche( Base ):
    functionname = 'rename_roche'

    def test_renames_roche_byregion( self ):
        files = self.mock_roche( 'R_2001_01_01_01_01_01_FLX00000001_adminrig_000000_TEST1', 'sample1', byregion=True )
        rw, rd, rbs = files
        self.print_tempdir()
        self._C( rbs, 'sample1', 'sample_1', 'NGSData' )
        self.print_tempdir()

        for f in files:
            ok_( not exists(f), '{} still exists'.format(f) )
        # Only need to do one exist here since it is a symlink all the way through
        newp = rbs.replace( 'sample1', 'sample_1' )
        ok_( exists(newp), '{} does not exist'.format(rbs) )

    def test_renames_roche_noregion( self ):
        files = self.mock_roche( 'R_2001_01_01_01_01_01_FLX00000001_adminrig_000000_TEST1', 'sample1', byregion=False )
        rw, rd, rbs = files
        self.print_tempdir()
        self._C( rbs, 'sample1', 'sample_1', 'NGSData' )
        self.print_tempdir()
        for f in files:
            ok_( not exists(f), '{} still exists'.format(f) )
        # Only need to do one exist here since it is a symlink all the way through
        newp = rbs.replace( 'sample1', 'sample_1' )
        ok_( exists(newp), '{} does not exist'.format(rbs) )

    def test_renames_roche_symlinkrawdata( self ):
        files = self.mock_roche( 'R_2001_01_01_01_01_01_FLX00000001_adminrig_000000_TEST1', 'sample1', byregion=False, rbsraw=True )
        rw, rd, rbs = files
        self.print_tempdir()
        self._C( rbs, 'sample1', 'sample_1', 'NGSData' )
        self.print_tempdir()
        for f in files:
            ok_( not exists(f), '{} still exists'.format(f) )
        # Only need to do one exist here since it is a symlink all the way through
        newp = rbs.replace( 'sample1', 'sample_1' )
        ok_( exists(newp), '{} does not exist'.format(rbs) )

class TestRenameMiSeq( Base ):
    functionname = 'rename_miseq'

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
    functionname = 'rename_sanger'

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

class TestRenameSample( Base ):
    functionname = 'rename_sample'

    def test_removes_empty_dirs( self ):
        self.mock_sanger( 'Run_3130xl_something', 'sample_1' )
        self.mock_roche( 'R_2013_10_11_16_42_14_FLX12070283_Administrator_10112013_Den1', 'sample_1' )
        self.mock_miseq( '010101_M00001_0001_000000000-A6F07', 'sample_1' )
        self.mock_iontorrent( '2001_01_01', 'sample_1' )

        self._C( 'sample_1', '00001', 'NGSData' )

        oldrbs = join( 'NGSData', 'ReadsBySample', 'sample_1' )
        newrbs = oldrbs.replace('sample_1', '00001')
        ok_( not exists(oldrbs), 'Did not remove empty directory sample_1' )
        ok_( exists(newrbs) and os.listdir(newrbs), 'Did not rename correctly to {}'.format(newrbs) )
