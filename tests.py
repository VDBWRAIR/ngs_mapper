from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import Mock, patch

import tempfile
import shutil
import os

reads = {
    'Roche454': [
        '1514A00101984N_CP2__1__TI5__2010_03_19__pH1N1.sff',
    ],
    'IonTorrent': [
        'Denv1_45AZ5_Virus_Lot1_82__1__IX031__2013_09_23__Den1.sff',
    ],
    'Sanger': [
        '1517010460_R1425_2011_08_24_H3N2_HA_0151_D09.fastq',
    ],
    'MiSeq': [
        '1090-01_S22_L001_R1_001_2013_12_10.fastq',
    ]
}

class BaseClass( object ):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp(prefix='test')
        os.chdir(self.tempdir)

    def tearDown(self):
        os.chdir('/')
        shutil.rmtree(self.tempdir)

    def mock_reads_dir(self,platforms=reads.keys()):
        readsdir = self.tempdir
        listing = {}
        for plat,readfiles in reads.items():
            listing[plat] = sorted(readfiles)
            # Create the files
            for readfile in listing[plat]:
                open(readfile,'w').close()
        return listing

    def reads_for_platforms(self):
        mapping = {}
        for plat, readfiles in reads.items():
            for readfile in readfiles:
                mapping[readfile] = plat
        return mapping

class TestFunctional(BaseClass):
    def test_reads_by_plat_emptydir(self):
        from data import reads_by_plat as rbp
        # no files so empty dict returned
        eq_( {}, rbp( self.tempdir ) )

    def test_reads_by_plat(self):
        import data
        with patch('data.platform_for_read', lambda filepath: self.reads_for_platforms()[filepath]) as a:
            with patch('data.filter_reads_by_platform', lambda path,plat: reads[plat]) as b:
                for plat, readfiles in reads.items():
                    eq_( readfiles, sorted(data.reads_by_plat( '' )[plat]) )

    def test_filter_reads_by_platform_all(self):
        from data import filter_reads_by_platform as frbp
        expected = self.mock_reads_dir()
        for plat, readfiles in expected.items():
            result = frbp( self.tempdir, plat )
            eq_( readfiles, sorted(result) )

    def test_platform_for_read( self ):
        from data import platform_for_read as pfr
        with patch('data.platform_for_read', lambda filepath: self.reads_for_platforms()[filepath]):
            for readfile, plat in self.reads_for_platforms().items():
                eq_( plat, pfr( readfile ) )

    def test_no_platform_for_read( self ):
        from data import platform_for_read as pfr
        from data import NoPlatformFound
        try:
            pfr( 'im_lonely' )
            assert False, "Did not raise NoPlatformFound"
        except NoPlatformFound as e:
            pass

class TestIntegration(BaseClass):
    def test_platform_for_read( self ):
        from data import platform_for_read as pfr
        for readfile, plat in self.reads_for_platforms().items():
            eq_( plat, pfr( readfile ) )

    def test_reads_by_plat(self):
        from data import reads_by_plat as rdp
        # Test each platform's file format
        expected = self.mock_reads_dir()
        result = rdp( self.tempdir )
        for plat,readfiles in expected.items():
            assert plat in result, "platform {} not in result".format(plat)
            eq_( readfiles, sorted(result[plat]) )
