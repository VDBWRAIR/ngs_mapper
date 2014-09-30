from imports import *

from miseqpipeline.data import NoPlatformFound

class Base(common.BaseClass):
    modulepath = 'miseqpipeline.data'

    def setUp(self):
        super(Base,self).setUp()
        self.reads = {
            'Roche454': [
                '1514A00101984N_CP2__1__TI5__2010_03_19__pH1N1.sff',
                'Vero_WHO_S16803__RL45__2012_05_01__Den2.fastq'
            ],
            'IonTorrent': [
                'Denv1_45AZ5_Virus_Lot1_82__1__IX031__2013_09_23__Den1.sff',
            ],
            'Sanger': [
                '1517010460_R1425_2011_08_24_H3N2_HA_0151_D09.fastq',
                '1517010460_F1425_2011_08_24_H3N2_HA_0151_D09.fastq',
            ],
            'MiSeq': [
                '1090-01_S22_L001_R1_001_2013_12_10.fastq',
                '1090-01_S22_L001_R2_001_2013_12_10.fastq',
            ]
        }

    def tearDown(self):
        super(Base,self).tearDown()

    def mock_reads_dir(self,path,platforms=None):
        if platforms is None:
            platforms = self.reads.keys()
        readsdir = self.tempdir
        listing = {}
        for plat,readfiles in self.reads.items():
            listing[plat] = sorted([join(path,r) for r in readfiles])
            # Create the files
            for readfile in listing[plat]:
                open(readfile,'w').close()
        return listing

    def reads_for_platforms(self):
        mapping = {}
        for plat, readfiles in self.reads.items():
            for readfile in readfiles:
                mapping[readfile] = plat
        return mapping

class TestIsSangerReadFile(Base):
    functionname = 'is_sanger_readfile'

    def test_detects_sanger( self ):
        f = 'sample1_F1_1979_01_01_Den2_Den2_0001_A01.fastq'
        f = join( fixtures.FIXDIR, 'trim_reads', f )
        r = self._C( f )
        ok_( r, 'Did not detect Sanger read' )

    def test_detects_not_sanger( self ):
        f = '1121__2__TI86__2012_04_04__Den2.fastq'
        f = join( fixtures.FIXDIR, 'trim_reads', f )
        r = self._C( f )
        ok_( not r, 'Incorrectly detected another platform as Sanger' )

class TestUnitPlatformForRead(Base):
    functionname = 'platform_for_read'

    def test_idents_platform( self ):
        for readfile, plat in self.reads_for_platforms().items():
            eq_( plat, self._C( readfile ) )

    def test_no_platform_for_read( self ):
        try:
            self._C( 'im_lonely' )
            assert False, "Did not raise NoPlatformFound"
        except NoPlatformFound as e:
            pass

    @raises(NoPlatformFound)
    def test_reads_by_plat_unkownfiles(self):
        self.reads['Sanger'].append( '1517010460_R1425_2011_08_24_H3N2_HA_0151_D09.junk' )
        self.test_idents_platform()

    def test_sff_now_fastq( self ):
        # Ensure that names that work with .sff exension also work with .fastq extension
        for plat, reads in self.reads.items():
            if plat in ('Roche454','IonTorrent'):
                for read in reads:
                    eq_( plat, self._C( read ) )
                    eq_( plat, self._C( read.replace('.sff','.fastq') ) )

class TestFunctional(Base):
    def test_reads_by_plat_emptydir(self):
        from miseqpipeline.data import reads_by_plat as rbp
        # no files so empty dict returned
        eq_( {}, rbp( self.tempdir ) )


    def test_reads_by_plat(self):
        with patch('miseqpipeline.data.platform_for_read', lambda filepath: self.reads_for_platforms()[filepath]) as a:
            with patch('miseqpipeline.data.filter_reads_by_platform', lambda path,plat: self.reads[plat]) as b:
                from miseqpipeline import data
                for plat, readfiles in self.reads.items():
                    if plat == 'MiSeq':
                        readfiles = [tuple(readfiles),]
                    eq_( sorted(readfiles), sorted(data.reads_by_plat( '' )[plat]) )

    def test_filter_reads_by_platform_all(self):
        from miseqpipeline.data import filter_reads_by_platform as frbp
        expected = self.mock_reads_dir(self.tempdir)
        for plat, readfiles in expected.items():
            result = frbp( self.tempdir, plat )
            eq_( readfiles, sorted(result) )

    def test_filter_reads_by_platform_bad_read_skipped(self):
        from miseqpipeline.data import filter_reads_by_platform as frbp
        fn = 'Sample_F1_1979_01_01_Den1_Den1_0001_A01.junk' 
        self.reads['Sanger'].append( fn )
        expected = self.mock_reads_dir(self.tempdir)
        for plat, readfiles in expected.items():
            result = frbp( self.tempdir, plat )
            if plat == 'Sanger':
                # Remove the read as it should not exist
                del readfiles[readfiles.index( join(self.tempdir,fn) )]
            eq_( readfiles, sorted(result) )

    def test_platform_has_paired_reads(self):
        reads = ['read1_r1','read1_r2','read2_r1','read3_r1','read3_r2']
        mates = [1,-1,4]
        with patch('miseqpipeline.data.find_mate', MagicMock(side_effect=mates)):
            from miseqpipeline.data import pair_reads
            expected = [('read1_r1','read1_r2'),'read2_r1',('read3_r1','read3_r2')]
            eq_( expected, pair_reads(reads) )

    def test_find_mate_read_has_mate(self):
        from miseqpipeline.data import find_mate
        eq_( 1, find_mate(self.reads['MiSeq'][0], self.reads['MiSeq']) )

    def test_find_mate_single_item_list(self):
        from miseqpipeline.data import find_mate
        eq_( -1, find_mate(self.reads['Roche454'][0], self.reads['Roche454']) )

    def test_find_mate_does_not_have_mate(self):
        from miseqpipeline.data import find_mate
        eq_( -1, find_mate(self.reads['Sanger'][0], self.reads['Sanger']) )

class TestIntegration(Base):
    def test_platform_for_read( self ):
        from miseqpipeline.data import platform_for_read as pfr
        for readfile, plat in self.reads_for_platforms().items():
            eq_( plat, pfr( readfile ) )

    def test_reads_by_plat_individual(self):
        from miseqpipeline.data import reads_by_plat as rdp
        # Test each platform's file format
        expected = self.mock_reads_dir(self.tempdir)
        result = rdp( self.tempdir )
        for plat,readfiles in expected.items():
            assert plat in result, "platform {} not in result".format(plat)
            if plat == 'MiSeq':
                readfiles = [tuple(sorted(readfiles)),]
                eq_( readfiles, [tuple(sorted(x)) for x in result[plat]] )
            else:
                eq_( readfiles, sorted(result[plat]) )

    def test_reads_by_all_unkownformats_are_skipped(self):
        self.reads['Sanger'].append( '1517010460_R1425_2011_08_24_H3N2_HA_0151_D09.junk' )
        expected = self.mock_reads_dir(self.tempdir)
        # remove the ab1 file from the expected as it should have been skipped
        expected['Sanger'] = [f for f in expected['Sanger'] if not f.endswith('.junk')]
        print expected['Sanger']
        self.mock_expected( expected )

    def test_reads_by_all(self):
        # Test each platform's file format
        expected = self.mock_reads_dir(self.tempdir)
        self.mock_expected( expected )

    def mock_expected( self, expected ):
        from miseqpipeline.data import reads_by_plat as rdp
        result = rdp( self.tempdir )
        ex = {}
        for plat, reads in expected.items():
            if plat == 'MiSeq':
                eq_( sorted(reads), sorted(result[plat][0]) )
            else:
                eq_( sorted(reads), sorted(result[plat])  )
        eq_( expected.keys(), result.keys() )

    def test_platform_has_paired_reads(self):
        from miseqpipeline.data import pair_reads
        eq_( [tuple(self.reads['MiSeq']),], pair_reads( self.reads['MiSeq'] ) )

    def test_platform_does_not_have_paired_reads(self):
        from miseqpipeline.data import pair_reads
        eq_( self.reads['Sanger'], pair_reads( self.reads['Sanger'] ) )
