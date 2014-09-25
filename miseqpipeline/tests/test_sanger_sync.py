from imports import *
import common
import string

class Base( common.BaseClass ):
    sampleabi = join( fixtures.FIXDIR, 'sample.ab1' )
    rund = 'Run_3130xl_{year}-{month}-{day}_01-01_{run}_{year}-{month}-{day}'
    sample = '{sample}_{primer}_{year}_{month}_{day}_{virus}_{gene}_{run}_{well}.{ext}'

    def _syncrun( self, ngsdata ):
        from sanger_sync import sync_run
        rund, samples, abilist = self.make_from_dir()
        r = sync_run( rund, ngsdata )
        rawd = join( ngsdata, 'RawData', 'Sanger', rund )
        return r, abilist, rund, rawd

    def _syncrunread( self, ngsdata ):
        from sanger_sync import sync_readdata
        r, abilist, rund, rawd = self._syncrun( ngsdata )
        sync_readdata( rawd, ngsdata )
        readdata = join( ngsdata, 'ReadData', 'Sanger', rund )
        return r, abilist, rund, rawd, readdata

    def make_from_dir( self, path=None ):
        if path is None:
            path = self.rund.format( year=2001, month=01, day=01, run=0001 )
        os.makedirs( path )

        samples = [
            {'sample':'sample1','primer':'F1','virus':'Den3','gene':'Den3','run':0001,'well':'A01','ext':'ab1'},
            {'sample':'sample_2','primer':'R1','virus':'Den3','gene':'Den3','run':0001,'well':'B01','ext':'ab1'},
            {'sample':'sample-3','primer':'F1','virus':'FluB','gene':'PB2','run':0001,'well':'C01','ext':'ab1'},
        ]

        abilist = []
        for sample in samples:
            sample['year'] = '2001'
            sample['month'] = '01'
            sample['day'] = '01'
            filename = self.sample.format( **sample )
            abilist.append( filename )
            filepath = join( path, filename )

            shutil.copy( self.sampleabi, filepath )

        return path, samples, abilist

    def check_rawdata( self, abilist, rawd ):
        ''' abilist is just basename list of abi files expected '''
        abi = set( os.listdir( rawd ) )
        abilist_set = set( abilist )
        eq_( set([]), abilist_set - abi, 'Abi lists were not identical {} != {}'.format(abi, abilist_set) )

    def check_readdata( self, abilist, readd ):
        ''' abilist is just basename list of abi files expected '''
        abi = [basename(r) for r in glob( join( readd, '*.ab1' ) )]
        abi = set( abi )
        abilist_set = set( abilist )
        eq_( set([]), abilist_set - abi, 'Abi lists were not identical {} != {}'.format(abi, abilist_set) )
        for abi in abilist:
            a = join( readd, abi )
            ok_( islink( a ), '{} is not a valid file/symlink'.format(a) )
            j = a.replace('.ab1', '.fastq')
            ok_( isfile( j ), 'Fastq version of {} does not exist'.format(j) )

    def print_ngs( self ):
        print subprocess.check_output( 'find . -ls', shell=True )

    def check_rbs( self, abilist, rbs ):
        ''' abilist is just basename list of abi files expected '''
        from sanger_sync import samplename_from_read
        self.print_ngs()
        for abi in abilist:
            sn = samplename_from_read( abi )
            rbsd = join( rbs, sn )
            read = join( rbsd, abi )
            ok_( islink( read ), '{}({}) is not a valid symlink'.format(read,os.readlink(read)) )
            fastq = read.replace('.ab1','.fastq')
            ok_( islink( fastq ), '{}({}) is not a valid symlink'.format(fastq,os.readlink(read)) )

class TestSync( Base ):
    def _C( self, *args, **kwargs ):
        from sanger_sync import sync_sanger
        return sync_sanger( *args, **kwargs )

    def test_sync_works( self ):
        rund, samples, abilist = self.make_from_dir()
        ngsdata = 'NGSData'
        self._C( rund, ngsdata )
        self.check_rawdata( abilist, join(ngsdata,'RawData','Sanger',rund) )
        self.check_rawdata( abilist, join(ngsdata,'ReadData','Sanger',rund) )
        self.check_rbs( abilist, join(ngsdata,'ReadsBySample') )

class TestFunctional( Base ):
    def _C( self, args ):
        script = TestFunctional.script_path( 'sanger_sync.py' )
        cmd = script + ' ' + args
        return TestFunctional.run_script( cmd )

    def test_runscorrectly( self ):
        rund, samples, abilist = self.make_from_dir()
        ngsdata = 'NGSData'
        args = '--ngsdata {} {}'.format(ngsdata, rund)
        retcode, output = self._C( args )
        eq_( 0, retcode, output )
        self.print_ngs()
        self.check_rawdata( abilist, join(ngsdata,'RawData','Sanger',rund) )
        self.check_rawdata( abilist, join(ngsdata,'ReadData','Sanger',rund) )
        self.check_rbs( abilist, join(ngsdata,'ReadsBySample') )

class TestSampleNameFromRead( Base ):
    def _C( self, *args, **kwargs ):
        from sanger_sync import samplename_from_read
        return samplename_from_read( *args, **kwargs )

    def test_valid_samplename( self ):
        validnames = [
            'sample1',
            'sample-1',
            'sample_1',
        ]
        for name in validnames:
            for p in 'FR':
                sn = self.sample.format(
                    sample=name, primer=p+'123',
                    year='2001', month='01', day='01',
                    virus='Vir1', gene='Gen1',
                    run='0001', well='A01', ext='fastq'
                )
                r = self._C( sn )
                eq_( name, r )

    def test_invalid_samplename( self ):
        from sanger_sync import InvalidFormat
        for name, p in (('sample 1', 'F123'), ('sample 1', 'F123')):
            sn = self.sample.format(
                sample=name, primer=p,
                year='2001', month='01', day='01',
                virus='Vir1', gene='Gen1',
                run='0001', well='A01', ext='fastq'
            )
            try:
                self._C( sn )
                ok_( False, 'Did not raise InvalidFormat exception' )
            except InvalidFormat as e:
                ok_( True )

class TestLinkReads( Base ):
    def _C( self, *args, **kwargs ):
        from sanger_sync import link_reads
        return link_reads( *args, **kwargs )

    def test_syncs_non_existing( self ):
        ngsdata = 'NGSData'
        reads, abilist, rund, rawd, readd = self._syncrunread( ngsdata )
        self._C( readd, ngsdata )
        rbs = join( 'NGSData', 'ReadsBySample' )
        self.check_rbs( abilist, rbs )

    def test_syncs_existing( self ):
        from sanger_sync import samplename_from_read
        ngsdata = 'NGSData'
        reads, abilist, rund, rawd, readd = self._syncrunread( ngsdata )
        rbs = join( 'NGSData', 'ReadsBySample' )
        # Now make one of the directories and ensure it was not recopied
        read = join( readd, abilist[0] )
        sn = samplename_from_read( read )
        abi = join( rbs, sn, basename( read ) )
        fq = abi.replace( '.ab1', '.fastq' )
        touch('read.ab1')
        touch('read.fq')
        os.makedirs( dirname(abi) )
        os.symlink('read.ab1',abi)
        os.symlink('read.fq',fq)
        self._C( readd, ngsdata )
        self.check_rbs( abilist, rbs )

        # Should not have changed the link
        eq_( 'read.ab1', os.readlink(abi) )
        eq_( 'read.fq', os.readlink(fq) )

class TestSyncRead( Base ):
    def _C( self, *args, **kwargs ):
        from sanger_sync import sync_readdata
        return sync_readdata( *args, **kwargs )

    def test_syncs_non_existing( self ):
        ngsdata = 'NGSData'
        reads, abilist, rund, rawd = self._syncrun( ngsdata )
        readd = join( ngsdata, 'ReadData', 'Sanger', rund )
        self._C( rawd, ngsdata )
        self.check_readdata( abilist, readd )

    def test_syncs_existing( self ):
        ngsdata = 'NGSData'
        reads, abilist, rund, rawd = self._syncrun( ngsdata )
        readd = join( ngsdata, 'ReadData', 'Sanger', rund )
        self._C( rawd, ngsdata )
        self.check_readdata( abilist, readd )
        # Remove file and try again
        abi = glob( join(readd,'*.ab1') )[0]
        fastq = glob( join(readd,'*.fastq') )[-1]
        os.unlink( abi )
        os.unlink( fastq )
        self._C( rawd, ngsdata )
        self.check_readdata( abilist, readd )

class TestSyncRun( Base ):
    def _C( self, *args, **kwargs ):
        from sanger_sync import sync_run
        return sync_run( *args, **kwargs )

    def test_invalid_filename( self ):
        from sanger_sync import InvalidFormat
        rund, samples, abilist = self.make_from_dir()
        seqfile = self.sample.format( 
            sample='sample 1', primer='F1', year=2001,
            month=01, day=01, virus='Den1',
            gene='Den1', run=0001, well='H01', ext='ab1'
        )
        seqf = join( rund, seqfile )
        touch( seqf )
        try:
            self._C( rund, 'NGSData' )
            ok_( False, 'Did not raise InvalidFormat' )
        except InvalidFormat as e:
            ok_( True )

    def test_syncs_run_no_dir( self ):
        rund, samples, abilist = self.make_from_dir()
        # Put a .seq file in there too just to make sure
        #  it doesn't transfer
        seqfile = self.sample.format( 
            sample='sample1', primer='F1', year=2001,
            month=01, day=01, virus='Den1',
            gene='Den1', run=0001, well='H01', ext='seq'
        )
        seqf = join( rund, seqfile )
        touch( seqf )
        self._C( rund, 'NGSData' )
        self.check_rawdata( abilist, join('NGSData','RawData','Sanger',rund) )

    def test_syncs_run_dir_exist_with_files( self ):
        rund, samples, abilist = self.make_from_dir()
        rawd = join('NGSData','RawData','Sanger',rund)
        os.makedirs( rawd )
        existread = self.sample.format( 
            sample='existing', primer='F1', year=2001,
            month=01, day=01, virus='Den1',
            gene='Den1', run=0001, well='H01', ext='ab1'
        )
        abilist.append( existread )
        touch( join(rawd,existread) )
        s = os.stat( join(rawd,existread) )
        self._C( rund, 'NGSData' )
        self.check_rawdata( abilist, rawd )
        # Existing file was not changed
        s2 = os.stat( join(rawd,existread) )
        eq_( s.st_size, s2.st_size )
        eq_( s.st_ctime, s2.st_ctime )
        eq_( s.st_atime, s2.st_atime )
        eq_( s.st_mtime, s2.st_mtime )

    def test_check_returnvalue( self ):
        rund, samples, abilist = self.make_from_dir()
        rawd = join( 'NGSData', 'RawData', 'Sanger', rund )
        r = self._C( rund, 'NGSData' )
        expect = set([join(rawd,abi) for abi in abilist])
        result = set(r)
        eq_( set([]), expect ^ result, 'Return value {} != {} Diff: {}'.format(result,expect,expect ^ result) )
