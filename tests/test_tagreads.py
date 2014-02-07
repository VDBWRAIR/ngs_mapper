from os.path import *
import os
import sys
from glob import glob
import subprocess
import tempfile
import shutil

from nose.tools import eq_, raises

class Base(object):
    def setUp( self ):
        self.tempdir = tempfile.mkdtemp( prefix='testTagReads' )
        this = dirname( __file__ )
        self.fixturedir = join( this, 'fixtures', 'tagreads' )
        self.bam = join( self.fixturedir, 'sample1.untagged.bam' )
        self.bai = join( self.fixturedir, 'sample1.untagged.bam.bai' )
        #sai = join( fixturedir, 'sample1.untagged.sai' )

        self.read_group_ids = ('Roche454','IonTorrent','Sanger','MiSeq')
        self.read_group_platforms = ('L454','IONTORRENT','CAPILLARY','ILLUMINA')

    def tearDown( self ):
        shutil.rmtree( self.tempdir )

    def copy_to_temp( self, filepath ):
        # Copy file to tempdir and return new path
        tmppath = join( self.tempdir, basename( filepath ) )
        shutil.copy( filepath, tmppath )
        return tmppath

    def temp_copy_files( self ):
        # Can call this so we don't copy every test by putting it in
        # setUp
        self.bam = self.copy_to_temp( self.bam )
        self.bai = self.copy_to_temp( self.bai )

    def get_readgroup( self, bam, rg ):
        cmd = ['samtools', 'view', '-r', rg, bam]
        return subprocess.check_output( cmd )

    def get_numreads( self, bam ):
        import pysam
        s = pysam.Samfile( bam )
        nreads = s.mapped + s.unmapped
        s.close()
        return nreads

    def check_tagreadcounts( self, bamfile ):
        read_count = 0
        for rg in self.read_group_ids:
            c = len( self.get_readgroup( bamfile, rg ).splitlines() )
            read_count += c
            print "Read Group: {} had {} reads".format(rg,c)

        eq_( self.get_numreads(bamfile), read_count )

class TestUnitGetRGForRead(Base):
    def _C( self, read ):
        from tagreads import get_rg_for_read
        return get_rg_for_read( read )

    def maread( self, **kwargs ):
        from mock import MagicMock
        aread = MagicMock( 'pysam.AlignedRead' )
        for k,v in kwargs.items():
            setattr( aread, k, v )
        return aread

    def test_roche( self ):
        rn = 'IPAOEAB02G8KRJ'
        mr = self.maread( qname=rn )
        rg = self._C( mr )
        eq_( 'Roche454', rg )

    def test_sanger( self ):
        rn = 'sample1_R2882_01_14_2014'
        mr = self.maread( qname=rn )
        rg = self._C( mr )
        eq_( 'Sanger', rg )

    def test_miseq( self ):
        rn = (
            'M02261:4:000000000-A6FWH:1:2106:7558:24138',
            'M02261:4:000000000-A6FWH:1:2106:75581:24138',
            'M02261:4:000000000-A6FWH:1:2106:7558:241381',
        )
        for r in rn:
            mr = self.maread( qname=r )
            rg = self._C( mr )
            eq_( 'MiSeq', rg )

    def test_iontorrent( self ):
        rn = 'E6I0Z:00005:00031'
        mr = self.maread( qname=rn )
        rg = self._C( mr )
        eq_( 'IonTorrent', rg )

class TestUnitGetHeader(Base):
    def _C( self, bamfile ):
        from tagreads import get_bam_header
        return get_bam_header( bamfile )

    def test_gets_header( self ):
        hdr = self._C( self.bam )
        assert 'SQ' in hdr, 'SQ key was not in bam header'
        assert 'HD' in hdr, 'HD key was not in bam header'

class TestUnitTagReads(Base):
    def _C( self, bamfile, hdr ):
        from tagreads import tag_reads
        return tag_reads( bamfile, hdr )

    def test_sets_header( self ):
        from tagreads import get_bam_header
        hdr = get_bam_header( self.bam )
        # Make a copy as we will edit the file
        self.temp_copy_files()
        rg = {
            'SM': 'sample1',
            'ID': 'testrg',
            'PL': 'ILLUMINA',
            'CN': 'SeqCenter'
        }
        hdr['RG'] = [rg]
        # Self.bam is now the copied version after calling temp_copy_files
        self._C( self.bam, hdr )
        new_hdr = get_bam_header( self.bam )
        assert 'RG' in new_hdr, "RG header did not get added to header"
        for k,v in rg.items():
            eq_( v, new_hdr['RG'][0][k] )

        assert self.get_numreads( self.bam ) > 0, "No reads in bamfile"

    def test_bam_has_rg_header( self ):
        from tagreads import get_bam_header, HeaderExists
        self.temp_copy_files()
        hdr = get_bam_header( self.bam )
        rg = {
            'SM': 'sample1',
            'ID': 'testrg',
            'PL': 'ILLUMINA',
            'CN': 'SeqCenter'
        }
        hdr['RG'] = [rg]
        self._C( self.bam, hdr )
        try:
            self._C( self.bam, hdr )
            assert False, "Did not raise HeaderExists"
        except HeaderExists as e:
            assert True

    def count_rg( self, bam ):
        ''' Count how many of each uniq read group id '''
        import pysam
        s = pysam.Samfile( bam )
        counts = {}
        for aread in s.fetch():
            tags = dict( aread.tags )
            id = tags['RG']
            if id not in counts:
                counts[id] = 0
            counts[id] += 1
        return counts

    def test_readsin_eq_readsout( self ):
        from tagreads import get_rg_headers
        self.temp_copy_files()
        hdr = get_rg_headers( self.bam )
        enreads = self.get_numreads( self.bam )
        self._C( self.bam, hdr )
        rnreads = self.get_numreads( self.bam )
        # The read counts should be the same
        eq_( enreads, rnreads )

    def test_tags_reads( self ):
        from tagreads import get_rg_headers
        self.temp_copy_files()
        hdr = get_rg_headers( self.bam )
        enreads = self.get_numreads( self.bam )
        self._C( self.bam, hdr )
        counts = self.count_rg( self.bam )
        eq_( 1, counts['IonTorrent'] )
        eq_( 1, counts['Roche454'] )
        eq_( 1, counts['Sanger'] )
        eq_( 996, counts['MiSeq'] )

class TestUnitGetRGHeaders(Base):
    def _C( self, bamfile, SM=None, CN=None ):
        from tagreads import get_rg_headers
        return get_rg_headers( bamfile, SM, CN )

    def test_adds_platform_read_groups( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam )
        r = dict( zip( self.read_group_ids, self.read_group_platforms ) )
        seen_id = []
        for rg in hdr['RG']:
            id = rg['ID']
            pl = rg['PL']
            eq_( r[id], pl )
            seen_id.append( id )
        eq_( sorted(self.read_group_ids), sorted(seen_id) )

    def test_sets_sm_to_filename( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam )
        for rg in hdr['RG']:
            eq_( rg['SM'], 'sample1.untagged' )

    def test_sets_sm( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam, SM='sample1' )
        for rg in hdr['RG']:
            eq_( rg['SM'], 'sample1' )

    def test_sets_cn( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam, CN='seqcenter' )
        for rg in hdr['RG']:
            eq_( rg['CN'], 'seqcenter' )

class TestIntegrate(Base):
    def _C( self, bamfiles, options=[] ):
        this = dirname( dirname( __file__ ) )
        cmd = [join( this, 'tagreads.py' )] + bamfiles + options
        subprocess.check_call( cmd )

    def test_tags_reads( self ):
        self.temp_copy_files()
        self._C( [self.bam] )
        self.check_tagreadcounts( self.bam )

    def test_does_multiple_bams( self ):
        import pysam
        self.temp_copy_files()
        bam2 = join( self.tempdir, 'sample2.bam' )
        bai2 = join( self.tempdir, 'sample2.bam.bai' )
        shutil.copy( self.bam, bam2 )
        shutil.copy( self.bai, bai2 )
        self._C( [self.bam, bam2] )
        for b in [self.bam, bam2]:
            self.check_tagreadcounts( b )
            n = basename(b).replace('.bam','')
            s = pysam.Samfile( b )
            assert all( rg['SM'] == n for rg in s.header['RG'] ) == True, "Did not set {} as SM for {}".format(n,b)

    def test_samplename_argument( self ):
        import pysam
        self.temp_copy_files()
        self._C( [self.bam], ['-SM', 'sample1'] )
        s = pysam.Samfile( self.bam )
        rgs = s.header['RG']
        # Ensure each read group contains the samplename set
        for rg in rgs:
            eq_( rg['SM'], 'sample1' )

    def test_seqencingcenter_argument( self ):
        import pysam
        self.temp_copy_files()
        self._C( [self.bam], ['-CN', 'seqcenter'] )
        s = pysam.Samfile( self.bam )
        rgs = s.header['RG']
        # Ensure each read group contains the samplename set
        for rg in rgs:
            eq_( rg['CN'], 'seqcenter' )
