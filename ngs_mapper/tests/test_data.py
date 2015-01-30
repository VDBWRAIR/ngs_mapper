from imports import *
import gzip

from ngs_mapper.data import NoPlatformFound

class Base(common.BaseClass):
    modulepath = 'ngs_mapper.data'

    def setUp(self):
        super(Base,self).setUp()
        self.reads = {
            'Roche454': [
                'Vero_WHO_S16803__RL45__2012_05_01__Den2.fastq'
            ],
            'IonTorrent': [
                'Denv1_45AZ5_Virus_Lot1_82__1__IX031__2013_09_23__Den1.fastq',
            ],
            'Sanger': [
                '1517010460_R1425_2011_08_24_H3N2_HA_0151_D09.fastq',
                '1517010460_F1425_2011_08_24_H3N2_HA_0151_D09.fastq',
            ],
            'MiSeq': [
                '1090-01_S22_L001_R1_001_2013_12_10.fastq',
                '1090-01_S22_L001_R2_001_2013_12_10.fastq',
                '1090-01_S22_L001_R1_001_2013_12_10.fastq.gz',
                '1090-01_S22_L001_R2_001_2013_12_10.fastq.gz',
                '1090-01_S22_L001_R1_001.fastq',
                '1090-01_S22_L001_R2_001.fastq',
                '1090-01_S22_L001_R1_001.fastq.gz',
                '1090-01_S22_L001_R2_001.fastq.gz',
            ]
        }

        self.read_ids = {
            'Roche454': 'AAAAAAAAAAAAAA',
            'IonTorrent': 'IIIII:1:1',
            'Sanger': 'foo-BAR',
            'MiSeq': '09aZ_:1000:000000000-A1B11:1000:0001:0001:0001'
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
                if readfile.endswith('.gz'):
                    fh = gzip.open(readfile, 'wb')
                else:
                    fh = open(readfile,'w')
                fh.write('@{0}\nATGC\n+\n!!!!\n'.format(
                    self.read_ids[plat]
                ))
                fh.close()
        return listing

    def reads_for_platforms(self):
        mapping = {}
        for plat, readfiles in self.reads.items():
            for readfile in readfiles:
                mapping[readfile] = plat
        return mapping

class TestPairReads(Base):
    functionname = 'pair_reads'

    @attr('current')
    def test_pairs_reads(self):
        readlist = self.reads['MiSeq']
        r = self._C( readlist )
        e = [(readlist[i], readlist[i+1]) for i in range(0, len(readlist), 2)]
        eq_( e, r )

    def test_pairs_reads_abs_paths(self):
        readlist = [join('/dev/shm/tdir',read) for read in self.reads['MiSeq']]
        r = self._C( readlist )
        e = [(readlist[i], readlist[i+1]) for i in range(0, len(readlist),2)]
        eq_( e, r )

class TestFindMate(Base):
    functionname = 'find_mate'

    def test_finds_mate_sorted(self):
        readlist = self.reads['MiSeq']
        print readlist

        for i, read in enumerate(readlist):
            r = self._C(read, readlist)
            # Odd index should be the one after
            if i % 2 == 0:
                expectindex = i + 1
            else:
                expectindex = i - 1

            eq_(expectindex, r)

    def test_finds_mate_reverse_sorted(self):
        readlist = self.reads['MiSeq'][::-1]
        r = self._C( readlist[0], readlist )
        eq_( 1, r )
        r = self._C( readlist[1], readlist )
        eq_( 0, r )

    def test_not_miseq_read(self):
        r = self._C( self.reads['Sanger'][0], self.reads['Sanger'] )
        eq_( -1, r )

    def test_read_miseq_with_R3(self):
        readlist = [self.reads['MiSeq'][0].replace('_R1_','_R3_')]
        r = self._C( readlist[0], readlist )
        eq_( -1, r )

    def test_miseq_with_no_mate(self):
        readlist = [self.reads['MiSeq'][0]]
        r = self._C( readlist[0], readlist )
        eq_( -1, r )

    def test_finds_mate_for_not_readsbysample_standard(self):
        read_filenames = [
            '1090-01_S22_L001_R1_001.fastq',
            '1090-01_S22_L001_R2_001.fastq',
        ]

    def test_finds_mate_for_gzip_extension(self):
        read_filenames = [
            '1090-01_S22_L001_R1_001.fastq.gz',
            '1090-01_S22_L001_R2_001.fastq.gz',
        ]

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

class TestFunctional(Base):
    def test_reads_by_plat_emptydir(self):
        from ngs_mapper.data import reads_by_plat as rbp
        # no files so empty dict returned
        eq_( {}, rbp( self.tempdir ) )

    @attr('current')
    def test_reads_by_plat(self):
        with patch('ngs_mapper.data.platform_for_read', lambda filepath: self.reads_for_platforms()[filepath]) as a:
            with patch('ngs_mapper.data.filter_reads_by_platform', lambda path,plat: self.reads[plat]) as b:
                from ngs_mapper import data
                for plat, readfiles in self.reads.items():
                    if plat == 'MiSeq':
                        readfiles = [(readfiles[i], readfiles[i+1]) for i in range(0, len(readfiles), 2)]
                    eq_( sorted(readfiles), sorted(data.reads_by_plat( '' )[plat]) )

    @attr('current')
    def test_filter_reads_by_platform_all(self):
        from ngs_mapper.data import filter_reads_by_platform as frbp
        expected = self.mock_reads_dir(self.tempdir)
        for plat, readfiles in expected.items():
            result = frbp(self.tempdir, plat)
            print plat
            print result
            eq_( readfiles, sorted(result) )

    def test_filter_reads_by_platform_bad_read_skipped(self):
        from ngs_mapper.data import filter_reads_by_platform as frbp
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
        with patch('ngs_mapper.data.find_mate', MagicMock(side_effect=mates)):
            from ngs_mapper.data import pair_reads
            expected = [('read1_r1','read1_r2'),'read2_r1',('read3_r1','read3_r2')]
            eq_( expected, pair_reads(reads) )

    def test_find_mate_read_has_mate(self):
        from ngs_mapper.data import find_mate
        eq_( 1, find_mate(self.reads['MiSeq'][0], self.reads['MiSeq']) )

    def test_find_mate_single_item_list(self):
        from ngs_mapper.data import find_mate
        eq_( -1, find_mate(self.reads['Roche454'][0], self.reads['Roche454']) )

    def test_find_mate_does_not_have_mate(self):
        from ngs_mapper.data import find_mate
        eq_( -1, find_mate(self.reads['Sanger'][0], self.reads['Sanger']) )

class TestIntegration(Base):
    def test_reads_by_plat_individual(self):
        from ngs_mapper.data import reads_by_plat as rdp
        # Test each platform's file format
        expected = self.mock_reads_dir(self.tempdir)
        result = rdp(self.tempdir)
        print result
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
        from ngs_mapper.data import reads_by_plat as rdp
        result = rdp( self.tempdir )
        ex = {}
        for plat, reads in expected.items():
            if plat == 'MiSeq':
                eq_( sorted(reads), sorted(result[plat][0]) )
            else:
                eq_( sorted(reads), sorted(result[plat])  )
        eq_( expected.keys(), result.keys() )

    @attr('current')
    def test_platform_has_paired_reads(self):
        from ngs_mapper.data import pair_reads
        readlist = self.reads['MiSeq']
        e = [(readlist[i], readlist[i+1]) for i in range(0, len(readlist), 2)]
        eq_( e, pair_reads( self.reads['MiSeq'] ) )

    def test_platform_does_not_have_paired_reads(self):
        from ngs_mapper.data import pair_reads
        eq_( self.reads['Sanger'], pair_reads( self.reads['Sanger'] ) )


# Better practices here

import mock
import unittest

from .. import data

@mock.patch('__builtin__.open', mock.MagicMock())
@mock.patch.object(data, 'gzip', mock.MagicMock())
@mock.patch.object(data, 'SeqIO')
class TestPlatformForRead(unittest.TestCase):
    def test_no_reads_in_file(self, mseqio):
        mseqio.parse.return_value = iter([])
        self.assertRaises(data.NoPlatformFound, data.platform_for_read, 'file.fastq')

    def test_identifies_roche_reads(self, mseqio):
        records = [
            mock.Mock(id='ABC1ABC12ABCA1', seq='ATGC'),
            mock.Mock(id='CBAABC12ABCDA1', seq='ATGC'),
        ]
        for record in records:
            mseqio.parse.return_value = iter([record])
            r = data.platform_for_read('/path/to/read.sff')
            self.assertEqual('Roche454', r)

    def test_identifies_ion_reads(self, mseqio):
        records = [
            mock.Mock(id='IIIII:00001:00001', seq='ATGC'),
            mock.Mock(id='AAAAA:1:1', seq='ATGC'),
            mock.Mock(id='BBBBB:10:1', seq='ATGC'),
            mock.Mock(id='BBBBB:10:10', seq='ATGC'),
        ]
        for record in records:
            mseqio.parse.return_value = iter([record])
            r = data.platform_for_read('/path/to/read.sff')
            self.assertEqual('IonTorrent', r)

    def test_identifies_sanger_reads(self, mseqio):
        records = [
            mock.Mock(id='SangerRead', seq='ATGC'),
            mock.Mock(id='aZ_09-', seq='ATGC')
        ]
        for record in records:
            mseqio.parse.return_value = iter([record])
            r = data.platform_for_read('/path/to/read.ab1')
            self.assertEqual('Sanger', r)
            self.assertEqual('abi', mseqio.parse.call_args[0][1])

    def test_identifies_miseq_reads(self, mseqio):
        #@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
        #@[a-zA-Z0-9_]+:[0-9]+:[a-zA-Z0-9-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+ [0-9]+:Y|N:[0-9]+:[ATGC]+
        records = [
            mock.Mock(id='aZ0_:1:aZ0-:1:1:1:1', seq='ATGC'),
            mock.Mock(id='_aZ_:1000:000000000-A1B11:1000:0001:0001:0001', seq='ATGC'),
        ]
        for record in records:
            mseqio.parse.return_value = iter([record])
            r = data.platform_for_read('/path/to/read.sff')
            self.assertEqual('MiSeq', r)

    def test_handles_gzipped_filepath(self, mseqio):
        record = mock.Mock(id='aZ_09-', seq='ATGC')
        mseqio.parse.return_value = iter([record])
        r = data.platform_for_read('/path/to/read.fastq.gz')
        self.assertEqual('Sanger', r)

    def test_raises_no_platform_found_when_no_match_is_found(self, mseqio):
        records = [
            mock.Mock(id='!#$%^&*(', seq='ATGC')
        ]
        for record in records:
            mseqio.parse.return_value = iter([record])
            self.assertRaises(
                data.NoPlatformFound,
                data.platform_for_read, '/path/to/read.sff'
            )

    def test_raises_no_platform_found_when_extension_not_supported_seqio(self, mseqio):
        mseqio.parse.side_effect = ValueError("Unknown format 'bar'")
        self.assertRaises(
            data.NoPlatformFound,
            data.platform_for_read, '/path/to/foo.bar'
        )

class TestFileHandle(unittest.TestCase):
    def test_returns_normal_handle(self):
        with mock.patch('__builtin__.open') as mopen:
            r = data.file_handle('/path/to/file.ext')
            self.assertEqual((mopen.return_value, 'ext'), r)

    def test_returns_normal_handle_basename(self):
        with mock.patch('__builtin__.open') as mopen:
            r = data.file_handle('file.ext')
            self.assertEqual((mopen.return_value, 'ext'), r)
        
    def test_returns_gzip_handle(self):
        with mock.patch.object(data, 'gzip') as mgzip:
            r = data.file_handle('/path/to/file.ext.gz')
            self.assertEqual((mgzip.open.return_value, 'ext'), r)
