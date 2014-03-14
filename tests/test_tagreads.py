from os.path import *
import os
import sys
from glob import glob
import subprocess
import tempfile
import shutil
import re

from nose.tools import eq_, raises, ok_
from nose.plugins.attrib import attr

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

    def mock_read_pysam( self ):
        from pysam import AlignedRead
        return AlignedRead()

    def mock_read( self ):
        from samtools import SamRow
        return SamRow('qname\t0\tRef1\t1\t60\t*\t=\t0\t0\t*\t*')

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
        read = self.mock_read()
        read.QNAME = rn
        rg = self._C( read )
        eq_( 'Roche454', rg )

    def test_sanger( self ):
        rn = 'sample1_R2882_01_14_2014'
        read = self.mock_read()
        read.QNAME = rn
        rg = self._C( read )
        eq_( 'Sanger', rg )

    def test_miseq( self ):
        rn = (
            'M02261:4:000000000-A6FWH:1:2106:7558:24138',
            'M02261:4:000000000-A6FWH:1:2106:75581:24138',
            'M02261:4:000000000-A6FWH:1:2106:7558:241381',
        )
        for r in rn:
            read = self.mock_read()
            read.QNAME = r
            rg = self._C( read )
            eq_( 'MiSeq', rg )

    def test_iontorrent( self ):
        rn = (
            'AAAAA:00000:00000',
            'AAAAA:00000:0',
            'AAAAA:0:00000',
            'AAAAA:0:0'
        )
        for r in rn:
            read = self.mock_read()
            read.QNAME = r
            rg = self._C( read )
            eq_( 'IonTorrent', rg )

class TestUnitGetHeader(Base):
    def _C( self, bamfile ):
        from tagreads import get_bam_header
        return get_bam_header( bamfile )

    def test_gets_header( self ):
        hdr = self._C( self.bam )
        ok_( isinstance( hdr, str ), 'Header is not a string' )
        ok_( 'SQ' in hdr, 'SQ key was not in bam header' )
        ok_( 'HD' in hdr, 'HD key was not in bam header' )

    def test_header_correct( self ):
        hdr = self._C( self.bam )
        with open( 't.sam', 'w' ) as fh:
            fh.write( hdr )
        import samtools
        h = samtools.view( 't.sam', S=True, H=True )
        eq_( hdr, h.read().rstrip() )

class TestUnitTagReads(Base):
    def _C( self, bamfile, hdr ):
        from tagreads import tag_reads
        return tag_reads( bamfile, hdr )

    def test_sets_header( self ):
        from tagreads import get_bam_header
        hdr = get_bam_header( self.bam )
        # Make a copy as we will edit the file
        self.temp_copy_files()
        # New read group header
        rg = '\n@RG\tID:testrg\tSM:sample1\tCN:SeqCenter\tPL:ILLUMINA\n'
        hdr += rg
        # Self.bam is now the copied version after calling temp_copy_files
        self._C( self.bam, hdr )
        new_hdr = get_bam_header( self.bam ).splitlines()
        reg = new_hdr[-1]
        eq_( rg.strip(), reg.strip(), 'Read group testrg did not make it into the header' )
        assert self.get_numreads( self.bam ) > 0, "No reads in bamfile"

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

class TestUnitTagReadGroup(Base):
    def _C( self, read ):
        from tagreads import tag_readgroup
        return tag_readgroup( read )

    def test_tags_sanger( self ):
        read = self.mock_read()
        read.QNAME = 'TestingAread_sanger'
        self._C( read )
        eq_( 'Sanger', read.TAGS[0][1] )
        eq_( 1, len(read.TAGS) )

    def test_tags_miseq( self ):
        read = self.mock_read()
        read.QNAME = 'M00000:0:000000000-AAAAA:0:0000:0000:0000'
        self._C( read )
        eq_( 'MiSeq', read.TAGS[0][1] )
        eq_( 1, len(read.TAGS) )

    def test_tags_roche( self ):
        read = self.mock_read()
        read.QNAME = 'A' * 14
        self._C( read )
        eq_( 'Roche454', read.TAGS[0][1] )
        eq_( 1, len(read.TAGS) )

    def test_tags_iontorrent( self ):
        read = self.mock_read()
        read.QNAME = 'AAAAA:00000:00000'
        self._C( read )
        eq_( 'IonTorrent', read.TAGS[0][1] )
        eq_( 1, len(read.TAGS) )

class TestUnitTagRead(Base):
    def _C( self, read, tags ):
        from tagreads import tag_read
        return tag_read( read, tags )

    def test_no_tags( self ):
        read = self.mock_read()
        tags = []
        r = self._C( read, tags )
        eq_( str(read), str(r) )

    def test_all_duplicate_tags( self ):
        # When all tags are duplicates
        read = self.mock_read()
        read._tags = 'RG:Z:Test'
        # Same tag
        tags = ['RG:Z:Test']
        r = self._C( read, tags )
        eq_( [('RG','Test')], read.TAGS )
        eq_( 'RG:Z:Test', read._tags )

    def test_skips_identical_tags( self ):
        # Don't append a tag that is identical to one that already exists
        read = self.mock_read()
        # Give mock read a tag
        read._tags = 'RG:Z:Test'
        # Now try to append tags, but with one of the tags being the tag
        # that already exists in the read to make sure it gets skipped, but the other
        # gets added to avoid the issue where tags could be duplicated many times if the user
        # runs the script 2x
        tags = ['RG:Z:Test','RG:Z:Test2']
        r = self._C( read, tags )
        # Should be the original tag and the other tag that was different
        eq_( [('RG','Test'),('RG','Test2')], read.TAGS )
        eq_( '\t'.join(tags), read._tags )
        eq_( r.TAGS, read.TAGS )

    def test_sets_tags( self ):
        read = self.mock_read()
        tags = ['RG:Z:Test','RX:Z:Test2']
        r = self._C( read, tags )
        eq_( [('RG','Test'),('RX','Test2')], read.TAGS )
        eq_( '\t'.join(tags), read._tags )
        eq_( r.TAGS, read.TAGS )

    def test_adds_tags( self ):
        read = self.mock_read()
        read._tags = 'RG:Z:Test1'
        tags = ['RG:Z:Test','RX:Z:Test2']
        r = self._C( read, tags )
        eq_( 'RG:Z:Test1\t' + '\t'.join(tags), read._tags )
        eq_( r.TAGS, read.TAGS )
        eq_( [('RG','Test1'),('RG','Test'),('RX','Test2')], read.TAGS )

    def test_skips_secondary_supplementary( self ):
        read = self.mock_read()
        read.FLAG = 2048
        tags = ['RG:Z:Test']
        r = self._C( read, tags )
        eq_( [], r.TAGS )

class TestUnitGetRGHeaders(Base):
    def _C( self, bamfile, SM=None, CN=None ):
        from tagreads import get_rg_headers
        return get_rg_headers( bamfile, SM, CN )

    def test_does_not_readd_headers( self ):
        # Make sure headers that exist are not duplicated
        from tagreads import get_bam_header
        import samtools
        self.temp_copy_files()
        # Create a new bam file that the header already has a read group in it
        # as well as a new header
        hdr = get_bam_header( self.bam )
        hdr += '\n' + '@RG\tID:Test\tCN:cn\tSM:sm\tPL:ILLUMINA\n'
        hdr += '@RG\tID:Roche454\tSM:312\tPL:L454\n'
        hdr += '@RG\tID:IonTorrent\tSM:312\tPL:IONTORRENT\n'
        hdr += '@RG\tID:MiSeq\tSM:312\tPL:ILLUMINA\n'
        hdr += '@RG\tID:Sanger\tSM:312\tPL:CAPILLARY\n'
        # Read the pipe which should be sam input and output bam
        s = samtools.view( self.bam )
        # Put in new header and sam output after it
        with open('t.sam','w') as fh:
            # New header
            fh.write( hdr )
            # Then reads
            fh.write(s.read())
        # Now convert that pipe to bam
        b = samtools.view( 't.sam', h=True, S=True, b=True )
        # Write the bam output
        with open( 'hasrg.bam', 'wb' ) as bamfh:
            bamfh.write( b.read() )
        # Close the file handles
        b.close()
        s.close()
        # Now we have a bamfile with an existing header that we can test
        r = self._C( 'hasrg.bam', 'sm', 'cn' )
        # Make sure that the new header made it in and that the MiSeq header was not duplicated
        header_lines = r.splitlines(True)
        read_groups = [rg for rg in header_lines if rg.startswith('@RG')]
        num_miseq = 0
        num_test = 0
        for rg in read_groups:
            if 'ID:Test\t' in rg:
                num_test += 1
            if 'ID:MiSeq\t' in rg:
                num_miseq += 1

        eq_( 1, num_miseq, "Header was duplicated which is incorrect" )
        eq_( 1, num_test, "Existing header was removed somehow" )
        # How many platform readgroups to expect(includes the MiSeq one that we are testing)
        i = len(self.read_group_ids)
        # Now increase that by HD, SQ and Test RG
        i += 3
        eq_( i, len(header_lines), "Incorrect number of header lines" )

    def test_correct_header( self ):
        # Make sure the expected number of header lines exist
        self.temp_copy_files()
        from tagreads import get_bam_header
        hdr = get_bam_header( self.bam )
        numlines = len(hdr.splitlines())
        hdr = self._C( self.bam, 'sample1', 'seqcenter' )
        linesexpected = len(self.read_group_ids) + numlines
        resultlines = len(hdr.splitlines())
        eq_( linesexpected, resultlines, "Expected {} header lines but got {}".format(linesexpected,resultlines) )

    def test_adds_platform_read_groups( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam )
        r = dict( zip( self.read_group_ids, self.read_group_platforms ) )
        seen_id = []
        hdr = hdr.splitlines()
        readgroups = [l for l in hdr if l.startswith('@RG')]
        for rg in readgroups:
            p = re.findall( '@RG\tID:(\S+)\tSM:(\S+)(?:\tCN:(\S+)){0,1}\tPL:(\S+)', rg )[0]
            id = p[0]
            pl = p[3]
            eq_( r[id], pl )
            seen_id.append( id )
        eq_( sorted(self.read_group_ids), sorted(seen_id) )

    def test_sets_sm_to_filename( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam )
        readgroups = [l for l in hdr.splitlines() if l.startswith('@RG')]
        for rg in readgroups:
            ok_( 'sample1.untagged' in rg, 'Samplename was not set to filename correctly' )

    def test_sets_sm( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam, SM='sample1' )
        readgroups = [l for l in hdr.splitlines() if l.startswith('@RG')]
        for rg in readgroups:
            ok_( 'sample1' in rg )

    def test_sets_cn( self ):
        self.temp_copy_files()
        hdr = self._C( self.bam, CN='seqcenter' )
        readgroups = [l for l in hdr.splitlines() if l.startswith('@RG')]
        for rg in readgroups:
            ok_( 'seqcenter' in rg, 'Samplename was not set to filename correctly' )

class TestIntegrate(Base):
    def _C( self, bamfiles, options=[] ):
        this = dirname( dirname( __file__ ) )
        cmd = [join( this, 'tagreads.py' )] + bamfiles + options
        subprocess.check_call( cmd )

    def test_does_not_duplicate( self ):
        from tagreads import get_bam_header
        # Make sure if you run tagreads 2x it doesn't duplicate tags and headers
        self.temp_copy_files()
        self._C( [self.bam] )
        size_before = os.stat( self.bam ).st_size
        self._C( [self.bam] )
        size_after = os.stat( self.bam ).st_size
        self.check_tagreadcounts( self.bam )
        eq_( size_before, size_after, "Running tagreads.py 2x produced a different size bam file" )

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
        import samtools
        self.temp_copy_files()
        self._C( [self.bam], ['-SM', 'sample1'] )
        s = samtools.view( self.bam, H=True )
        rgs = s.readlines()
        s.close()
        # Ensure each read group contains the samplename set
        for rg in [r for r in rgs if r.startswith('@RG')]:
            ok_( 'SM:sample1\t' in rg, "Samplename did not make it into the headers" )

    def test_seqencingcenter_argument( self ):
        import samtools
        self.temp_copy_files()
        self._C( [self.bam], ['-CN', 'seqcenter'] )
        s = samtools.view( self.bam, H=True )
        rgs = s.readlines()
        s.close()
        # Ensure each read group contains the samplename set
        for rg in [r for r in rgs if r.startswith('@RG')]:
            ok_( 'CN:seqcenter\t' in rg, "Sequencing center did not make it into the headers" )

    def test_no_extra_files( self ):
        self.temp_copy_files()
        efiles = glob( '*' ) + glob( join( self.tempdir, '*' ) )
        self._C( [self.bam] )
        rfiles = glob( '*' ) + glob( join( self.tempdir, '*' ) )
        eq_( efiles, rfiles )
