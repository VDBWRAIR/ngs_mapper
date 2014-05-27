from imports import *
import common

class Base( common.BaseClass ):
    sampleabi = join( fixtures.FIXDIR, 'sample.ab1' )
    rund = 'Run_3130xl_{year}-{month}-{day}_01-01_{run}_{year}-{month}-{day}'
    sample = '{sample}_{primer}_{year}_{month}_{day}_{virus}_{gene}_{run}_{well}.{ext}'

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
            sample['year'] = 2001
            sample['month'] = 01
            sample['day'] = 01
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
        abi = [basename(r) for r in glob.glob( join( readd, '*.ab1' ) )]
        abi = set( abi )
        abilist_set = set( abilist )
        eq_( set([]), abilist_set - abi, 'Abi lists were not identical {} != {}'.format(abi, abilist_set) )
        for abi in abilist:
            ok( islink( abi ), '{} is not a valid file/symlink'.format(abi) )
            j = join( readd, abi.replace('.ab1', '.fastq') )
            ok_( islink( j ), 'Fastq version of {} does not exist'.format(abi) )

    def check_rbs( self, abilist, rbs ):
        ''' abilist is just basename list of abi files expected '''
        for abi in abilist:
            sn = re.search( '(\S+?)_[RF]\d', abi ).groups(1)[0]
            rbsd = join( rbs, sn )
            read = join( rbsd, abi )
            ok_( islink( read ), '{} is not a valid symlink'.format(read) )
            fastq = read.replace('.ab1','.fastq')
            ok_( islink( fastq ), '{} is not a valid symlink'.format(fastq) )

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

class TestSyncRead( Base ):
    def _C( self, *args, **kwargs ):
        from sanger_sync import sync_readdata
        return sync_readdata( *args, **kwargs )

class TestSyncRun( Base ):
    def _C( self, *args, **kwargs ):
        from sanger_sync import sync_run
        return sync_run( *args, **kwargs )

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
