import fixtures
import common

from nose.tools import eq_, ok_, raises
from nose.plugins.attrib import attr
from mock import Mock, MagicMock, patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

from base_caller import VCF_HEAD
import vcf
from vcf.model import _Record

from os.path import *
import sys
import os
from StringIO import StringIO

class VCFBase(common.BaseClass):
    def setUp( self ):
        super(VCFBase, self).setUp()
        self.writer,self.fp = self.vcf_file()
        
    def vcf_file( self, filepath=None ):
        if filepath is None:
            filepath = join( self.tempdir, 'test.vcf' )
        template = StringIO( VCF_HEAD.format('test.vcf') )
        template.name = basename(filepath)
        return vcf.Writer( open(filepath,'w'), template=vcf.Reader(template) ), filepath

    def vcf_record( self, *args, **kwargs ):
        return _Record( *args, **kwargs )

    def blank_record( self, *args, **kwargs ):
        from base_caller import blank_vcf_row
        return blank_vcf_row( *args, **kwargs )

    # quicker alias
    vr = vcf_record
    br = blank_record

class TestUnitIterRefs(VCFBase):
    def _C( self, vcffile ):
        from vcf_consensus import iter_refs
        return iter_refs( vcffile )

    def test_is_iterable_singleref( self ):
        # Just one entry
        self.writer.write_record( self.br( 'Ref1', 'A', 1 ) )
        self.writer.close()
        r = self._C( self.fp ).next()
        ok_( isinstance( r, SeqRecord ), "Did not return a Bio.SeqRecord instance" )

    def test_rows_are_correct_type( self ):
        # Mock Reference sequences
        refseq = 'ACGT'*10
        # Multiple references with single bases each
        for i in range(5):
            self.writer.write_record( self.br( 'Ref{}'.format(i), refseq, 1 ) )
        self.writer.close()
        print open(self.fp).read()
        for row in self._C( self.fp ):
            ok_( isinstance( row, SeqRecord ), "Did not return a Bio.SeqRecord instance"  )
            ok_( row.id.startswith( 'Ref' ), "Sequence id did not start with Ref" )

    @raises(StopIteration)
    def test_empty_vcf( self ):
        r = self._C( self.fp ).next()

    def test_correct_consensus( self ):
        ref = 'ACGTN*-'*10
        # Vcf will have 2 references Ref1 & Ref2
        for i in range( 1, 3 ):
            # Generates a vcf with one ref and all CB are equal to Ref base
            for p in range( len(ref) ):
                rec = self.br( 'Ref{}'.format(i), ref, p+1, ref[p] )
                self.writer.write_record( rec )
        self.writer.close()
        for i, row in enumerate( self._C( self.fp ) ):
            # Correct id field
            eq_( 'Ref{}'.format(i+1), row.id )
            # Correct sequence field
            eq_( ref, str(row.seq) )

class TestUnitWriteFasta(VCFBase):
    def _C( self, records, outputfile ):
        from vcf_consensus import write_fasta
        return write_fasta( records, outputfile )

    def mock_records( self, num ):
        ''' Just get some blank records to use '''
        # 6 records with ids 1 - 6
        for i in range( num ):
            s = SeqRecord(
                Seq( str(generic_dna.letters), generic_dna ),
                id='Rec{}'.format(i)
            )
            yield s

    def compare_seqrecords( self, expectedrecords, fastafile ):
        ''' Compare expected records with SeqRecords in fastafile '''
        seqs = SeqIO.index( fastafile, 'fasta' )
        # Make sure lengths are the same
        count = 0
        for eseq in expectedrecords:
            ok_( seqs.get( eseq.id, False ), "{} did not make it into consensus fasta file" )
            eq_( str(eseq.seq), str(seqs[eseq.id].seq), "Sequences did not match" )
            count += 1
        le = len(expectedrecords)
        ls = len(seqs)
        eq_( le, ls, "Did not have the same amount of records "\
                "in file({}) as in expected list({})".format(le,ls) )

        return len(seqs)
    csr = compare_seqrecords

    def test_does_not_write_empty( self ):
        def emptygen():
            yield
        of = join( self.tempdir, 'out.fasta' )
        r = self._C( emptygen(), of )
        eq_( 0, r )
        ok_( not exists( of ), "Created output fasta even though no records were provided" )

    def test_writes_single_record( self ):
        of = join( self.tempdir, 'out.fasta' )
        recs = self.mock_records( 1 )
        r = self._C( recs, of )
        eq_( 1, r )
        self.csr( list(self.mock_records( 1 )), of )

    def test_writes_multiple_record( self ):
        of = join( self.tempdir, 'out.fasta' )
        recs = self.mock_records( 6 )
        r = self._C( recs, of )
        eq_( 6, r )
        self.csr( list(self.mock_records(6)), of )

class TestUnitOutputDiff(VCFBase):
    def _C( self, vcfrecord ):
        from vcf_consensus import output_diff
        return output_diff( vcfrecord )

class TestIntegrate(VCFBase):
    def _C( self, vcffile, **kwargs ):
        import subprocess
        d = join( dirname( dirname( __file__ ) ) )
        script = join( d, 'vcf_consensus.py' )
        cmd = [script, vcffile]
        if kwargs.get( 'showdiffs', False ):
            cmd += ['--show-diffs']
        if kwargs.get( 'o', False ):
            cmd += ['-o', kwargs.get( 'o' )]

        try:
            return subprocess.check_call( cmd )
        except subprocess.CalledProcessError as e:
            print e
            print e.output
            raise e

    def test_outputname_not_provided( self ):
        self.writer.close()
        r = self._C( self.fp )
        ok_( exists( self.fp.replace('.vcf','.fasta') ), "Did not create correct output name" )

    def test_outputname_provided( self ):
        self.writer.close()
        r = self._C( self.fp, o=self.fp+'.fasta' )
        ok_( exists( self.fp+'.fasta' ), "Did not create correct output name" )
