from imports import *

class Base(common.BaseClass):
    def check_trimmed_against_record( self, trimmedrec, record, trim=True ):
        ''' Check trimmed seq against initial record '''
        # Useful values to check
        origseq = record.seq._data
        trimseq = trimmedrec.seq._data
        trimqual = trimmedrec._per_letter_annotations['phred_quality']
        origqual = record._per_letter_annotations['phred_quality']
        # Where does trimseq start in origseq
        start = origseq.index( trimseq )
        end = start + len(trimseq)
        cql = record.annotations['clip_qual_left']
        cqr = record.annotations['clip_qual_right']
        print 'Orig:' + origseq
        print 'Trim:' + trimseq
        if trim:
            eq_( start, cql, 'clip_qual_left {} and beginning of trim {} did not match'.format(cql,start) )
            eq_( end, cqr, 'clip_qual_right {} and end of trim {} did not match'.format(cqr,end) )
            eq_( cqr - cql, len(trimmedrec._per_letter_annotations['phred_quality']), 'Did not trim qualities' )
        else:
            eq_( origseq, trimseq, 'Sequence was trimmed when it should not have been' )
            eq_( origqual, trimqual, 'Qualities were trimmed when they should not have been' )

def rand_sff_seqs( ):
    ''' To aid in mocking Bio.SeqIO.parse '''
    for i in range( 10 ):
        yield common.rand_seqrec( 100, 0, 0, 5, 10 )

@patch( 'reads.SeqIO.parse' )
class TestSffsToFastq(Base):
    def setUp( self ):
        super( TestSffsToFastq, self ).setUp()

    def _C( self, *args, **kwargs ):
        from reads import sffs_to_fastq
        return sffs_to_fastq( *args, **kwargs )

    def _pre( self, bioparse ):
        self.f1records = [r for r in rand_sff_seqs()]
        self.f2records = [r for r in rand_sff_seqs()]
        self.records = self.f1records + self.f2records
        # parse should now return f1records then f2records
        bioparse.side_effect = [self.f1records,self.f2records]

    def _post( self, outfile, count, trim=True ):
        from Bio import SeqIO
        ok_( exists(outfile), 'Did not create out.fastq' )
        eq_( 20, count )
        newrecs = SeqIO.index(outfile,'fastq')
        for origr in self.records:
            self.check_trimmed_against_record( newrecs[origr.id], origr, trim )

    def test_trims( self, bioparse ):
        self._pre( bioparse )
        outfile = 'out.fastq'
        count = self._C( ['file1.sff','file2.sff'], outfile, True )
        self._post( outfile, count )

    def test_no_trims( self, bioparse ):
        self._pre( bioparse )
        outfile = 'out.fastq'
        count = self._C( ['file1.sff','file2.sff'], outfile, False )
        self._post( outfile, count, False )

class TestClipSeqRecord(Base):
    def setUp( self ):
        super(TestClipSeqRecord,self).setUp()
        self.seqlen = 100
        self.cal = 0
        self.car = 0
        self.cql = 5
        self.cqr = 10
        self.rec = common.rand_seqrec( 100, 0, 0, 5, 10 )

    def _C( self, *args, **kwargs ):
        from reads import clip_seq_record
        return clip_seq_record( *args, **kwargs )

    def test_trims( self ):
        from copy import copy
        trimrec = self._C( self.rec )
        self.check_trimmed_against_record( trimrec, self.rec )

    def test_bug_767( self ):
        seq = 'tcagtccaagctgcgatCGCCGTTTCCCAGTAGGTCTCGAGAGGGCTGGGGATAATCCCTTCTGGTGTGTTT'
        rec = common.rand_seqrec( seq, 0, 0, 17, 0 )
        trimrec = self._C( rec )
        eq_( 'CGCCGTTTCCCAGTAGGTCTCGAGAGGGCTGGGGATAATCCCTTCTGGTGTGTTT', trimrec.seq._data, 'Did not properly handle IonTorrent clipping' )

    def test_sanity_check_1( self ):
        seq = 'tcagctcgcgtgtcTACGGTAGCAGAGACTTGGTCTCTCTGATGGCTGGGTTGGTATCTTAttggcgctgacatgagtttgtacgtcgtcagaattgccatcctgtttctttcgaatttagccatatattggtggtcgtaagactcctaaagccgtagtcttcaaccgtccaacgagttccaagcttctgtttgtgttggtgtagaacccatgtcgtcagtgttcacctaatcgactgatggcgcaagggagcgattacgnnnnnnn'
        rec = common.rand_seqrec( seq, 0, 0, 14, 61 )
        trimrec = self._C( rec )
        eq_( 'TACGGTAGCAGAGACTTGGTCTCTCTGATGGCTGGGTTGGTATCTTA', trimrec.seq._data, 'Did not trim sequence as expected' )

@patch('bwa.seqio.concat_files')
class TestFunctionalCompileReads(Base):
    def test_reads_abspath_basepath_relpath(self,mock):
        from reads import compile_reads
        outputdir = 'output'
        outputdir = join(self.tempdir,outputdir)
        # Should return/create 3 files
        reads = [('p1_1.fastq','p1_2.fastq'),'/path/to/np1.fastq',('/path/to/p2_1.fastq.fastq','../p2_1.fastq.fastq'),'../np2.fastq']
        files = ['F.fq','R.fq','NP.fq']
        expected = {
            'F':join(outputdir,files[0]),
            'R':join(outputdir,files[1]),
            'NP':join(outputdir,files[2])
        }
        eq_( expected, compile_reads( reads, outputdir ) )
        eq_( len(mock.call_args_list), 3 )

    def test_compile_reads_paired_and_unpaired(self,mock):
        from reads import compile_reads
        outputdir = 'output'
        outputdir = join(self.tempdir,outputdir)
        # Should return/create 3 files
        reads = [('p1_1.fastq','p1_2.fastq'),'np1.fastq',('p2_1.fastq.fastq','p2_1.fastq.fastq'),'np2.fastq']
        files = ['F.fq','R.fq','NP.fq']
        expected = {
            'F':join(outputdir,files[0]),
            'R':join(outputdir,files[1]),
            'NP':join(outputdir,files[2])
        }
        eq_( expected, compile_reads( reads, outputdir ) )
        eq_( len(mock.call_args_list), 3 )

    def test_compile_reads_paired_only_single(self,mock):
        from reads import compile_reads
        outputdir = 'output'
        outputdir = join(self.tempdir,outputdir)
        # Should return/create 3 files
        reads = [('p1_1.fastq','p1_2.fastq')]
        files = ['F.fq','R.fq','NP.fq']
        expected = {
            'F':join(outputdir,files[0]),
            'R':join(outputdir,files[1]),
            'NP':None
        }
        eq_( expected, compile_reads( reads, outputdir ) )
        eq_( len(mock.call_args_list), 2 )

    def test_compile_reads_paired_only_multiple(self,mock):
        from reads import compile_reads
        outputdir = 'output'
        outputdir = join(self.tempdir,outputdir)
        # Should return/create 3 files
        reads = [('p1_1.fastq','p1_2.fastq'),('p2_1.fastq','p2_2.fastq')]
        files = ['F.fq','R.fq','NP.fq']
        expected = {
            'F':join(outputdir,files[0]),
            'R':join(outputdir,files[1]),
            'NP':None
        }
        eq_( expected, compile_reads( reads, outputdir ) )
        print mock.call_args_list
        eq_( len(mock.call_args_list), 2 )

    def test_compile_reads_unpaired_only_single(self,mock):
        from reads import compile_reads
        outputdir = 'output'
        outputdir = join(self.tempdir,outputdir)
        # Should return/create 3 files
        reads = ['p1.fastq']
        files = ['F.fq','R.fq','NP.fq']
        expected = {
            'F':None,
            'R':None,
            'NP':join(outputdir,files[2])
        }
        eq_( expected, compile_reads( reads, outputdir ) )
        eq_( len(mock.call_args_list), 1 )

    def test_compile_reads_unpaired_only_multiple(self,mock):
        from reads import compile_reads
        outputdir = 'output'
        outputdir = join(self.tempdir,outputdir)
        reads = ['p1.fastq','p2.fastq']
        files = ['F.fq','R.fq','NP.fq']
        expected = {
            'F':None,
            'R':None,
            'NP':join(outputdir,files[2])
        }
        eq_( expected, compile_reads( reads, outputdir ) )
        eq_( len(mock.call_args_list), 1 )

    @patch('reads.sffs_to_fastq')
    def test_converts_sff_to_fastq(self,sffs_to_fastq, concat_files):
        from reads import compile_reads
        def side_effect( ins, out ):
            fh = open(out,'w')
            fh.write('\n')
            fh.close()
        sffs_to_fastq.side_effect = side_effect
        outputdir = join(self.tempdir,'output')
        reads = ['np.sff','np.fastq',('F.fastq','R.fastq')]
        compile_reads( reads, outputdir )
        eq_( ['np.sff'], sffs_to_fastq.call_args_list[0][0][0] )
        eq_( 3, len(concat_files.call_args_list) )

    def test_compile_reads_non_supported_read_types(self,mock):
        # Only support fastq and sff right now
        from reads import compile_reads, InvalidReadFile
        outputdir = join(self.tempdir,'output')
        reads = ['np.ab1','np.fastq.gz']
        try:
            compile_reads( reads, outputdir )
            assert False, "Did not raise InvalidReadFile"
        except InvalidReadFile as e:
            pass
        except Exception as e:
            assert False, "Did not raise InvalidReadFile"

    def test_compile_reads_three_item_tuple(self,mock):
        from reads import compile_reads, InvalidReadFile
        outputdir = join(self.tempdir,'output')
        reads = [('one.fastq','two.fastq','three.fastq')]
        try:
            compile_reads( reads, outputdir )
            assert False, "ValueError not raised"
        except ValueError as e:
            pass

    def test_compile_reads_emptyreadfilelist(self,mock):
        from reads import compile_reads, InvalidReadFile
        outputdir = join(self.tempdir,'output')
        reads = ['np.sff','np.ab1','np.fastq.gz']
        eq_({'F':None,'R':None,'NP':None}, compile_reads( [], outputdir ) )

class TestUnitIsValidRead(object):
    def _C( self, readpath ):
        from reads import is_valid_read
        return is_valid_read( readpath )

    def test_invalid_read( self ):
        for ext in ('ab1','gz','adasdf'):
            eq_( False, self._C( 'file.'+ext ) )

    def test_valid_read( self ):
        from reads import VALID_READ_EXT
        for ext in VALID_READ_EXT:
            eq_( True, self._C( 'file.'+ext ) )

    def test_file_has_no_ext( self ):
        eq_( False, self._C( 'file' ) )

class TestIntegrationCompileReads(Base):
    def test_compile_reads_outputdir_not_exist(self):
        outputdir = join(self.tempdir,'output')
        self.run_compile_reads(outputdir)

    def test_compile_reads_outputdir_exists(self):
        outputdir = join(self.tempdir,'output')
        os.mkdir(outputdir)
        self.run_compile_reads(outputdir)

    def run_compile_reads(self,outputdir):
        from reads import compile_reads
        readsdir = join( fixtures.THIS, 'fixtures', 'reads' )
        sff = glob( join(readsdir,'*.sff') )
        miseq = tuple( glob( join(readsdir,'*L001*') ) )
        sanger = glob( join( readsdir, '*0001*' ) )
        reads = sff + [miseq] + sanger

        # Read dir will allow relative and abspath testing
        readdir = 'reads'
        os.mkdir(readdir)
        expected_linecounts = {}
        # Only 1 read F & R since just copied sanger read to fake the miseq paired
        expected_linecounts['F'] = 1*4
        expected_linecounts['R'] = 1*4
        # 100 sff reads + 2 sanger
        expected_linecounts['NP'] = 102*4

        result = compile_reads( reads, outputdir )
        assert os.access(result['F'],os.R_OK), "Did not create Forward file"
        assert os.access(result['R'],os.R_OK), "Did not create Reverse file"
        assert os.access(result['NP'],os.R_OK), "Did not create NonPaired file"
        for k,v in expected_linecounts.items():
            with open(result[k]) as fh:
                contents = fh.read()
                eq_( v, len(contents.splitlines()) )
