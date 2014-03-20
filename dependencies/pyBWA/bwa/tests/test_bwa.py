from nose.tools import eq_, raises
from ..bwa import BWA, BWAMem, BWAIndex
from .. import bwa
from .. import seqio

import tempfile
import shutil
import os
import os.path
import sys
import glob

import util
from util import get_bwa_path, INPUT_PATH, REF_PATH, BWA_PATH, test_bwa_available, ungzip

class BWASubTest( BWA ):
    ''' Test subclass to get around required_args exception '''
    def required_args( self ):
        pass

class BaseBWA( util.Base ):
    @classmethod
    def setUpClass( self ):
        super( BaseBWA, self ).setUpClass()
        self.bwa_path = self.mkbwa( 'test' )

    @classmethod
    def mkbwa( self, echomsg ):
        ''' Make fake bwa '''
        with open( 'bwa', 'w' ) as fh:
            fh.write( '#!/usr/bin/env bash\necho {}'.format(echomsg) )
        os.chmod( 'bwa', 0700 )
        return os.path.abspath( 'bwa' )

class TestBWA( BaseBWA ):
    @classmethod
    def setUpClass( self ):
        super( TestBWA, self ).setUpClass()
        self.rops = {rop:'value'+str(v) for v, rop in enumerate( BWA.REQUIRED_OPTIONS )}
        self.inst = BWASubTest( **self.rops )

    @raises( NotImplementedError )
    def test_requiredargs( self ):
        ''' Just make sure exception is raised '''
        BWA( **self.rops )

    def test_requireoptions_present( self ):
        ''' All required options present should pass '''
        BWASubTest( **self.rops )

    @raises( ValueError )
    def test_requireoptions_missing( self ):
        ''' Missing some required options should fail '''
        bwa = BWASubTest()

    def test_requiredoptions( self ):
        ''' Make sure required options are set as class vars '''
        bwa = self.inst
        erop = dict( zip( BWA.REQUIRED_OPTIONS, bwa.required_options_values ) )
        print bwa.required_options_values
        eq_( self.rops, erop )

    def test_bwareturncode_success( self ):
        ''' Just check to make sure not having usage statements works '''
        stderr = '''This does not have Ussage:  bbwa in it'''
        eq_( 0, self.inst.bwa_return_code( stderr ) )
        stderr = '''This only has Usage:  in it'''
        eq_( 0, self.inst.bwa_return_code( stderr ) )

    def test_bwareturncode_usage( self ):
        ''' Having Usage in the output should return failure '''
        stderr = '''This has Usage:bwa in it'''
        brc = self.inst.bwa_return_code
        eq_( 2, brc( stderr ) )
        stderr = '''This has Usage: bwa'''
        eq_( 2, brc( stderr ) )
        stderr = '''
this has Usage:  bwa
'''
        eq_( 2, brc( stderr ) )

    def test_runbwa_defaultoutput( self ):
        ''' Make sure default output works '''
        output_text = 'test output'
        bwa_path = self.mkbwa( '{}'.format(output_text) )
        # Get default filename argument from run_bwa function
        default_output_fn = self.inst.run_bwa.im_func.func_defaults[0]
        self.inst.run_bwa( [bwa_path], [], [] )
        eq_( os.path.exists( default_output_fn ), True )
        with open( default_output_fn ) as fh:
            eq_( output_text, fh.read().strip() )
        os.unlink( default_output_fn )

    def test_runbwa_optionsargs( self ):
        ''' Make sure options and args sent in are actually used '''
        bwa_path = self.mkbwa( '$@' )
        self.inst.run_bwa( [bwa_path], ['a','b'], ['c','d'], 'output' )
        with open( 'output' ) as fh:
            result_output = fh.read().strip()
        eq_( 'a b c d', result_output )

    def test_compilebwaoptions( self ):
        ''' Make sure options supplied to constructor make it through to running '''
        bwa_path = self.mkbwa( '$@' )
        ops = {'R':None, 'E':'', 't':3, 'a':5, 'f':True, 'j':1, 'bwa_path':bwa_path, 'command':'cmd'}
        bwa = BWASubTest( **ops )
        del ops['bwa_path']
        del ops['command']
        bwa.run( 'output' )
        with open( 'output' ) as fh:
            result_output = fh.read().strip()
        print "Expect: {}".format(ops)
        print "Result: {}".format(result_output)
        for op,val in ops.items():
            uval = ''
            val=str(val)
            if val.lower() not in ('true','false'):
                uval = ' ' + val
            if val.lower() not in ('','none'):
                opval = '-{}{}'.format(op,uval) 
                assert opval in result_output, '{} not in {}'.format(opval, result_output)
            else:
                assert op not in result_output, '{} should not be in {}'.format(op, result_output)

    @raises( ValueError )
    def test_runbwa_invalidpath( self ):
        ''' make sure path to bwa is valid '''
        self.inst.run_bwa( ['bad/path/to/bwa'], [], [] )

    def test_runbwa_invalidoptions( self ):
        ''' Integration test with bwa_return_code. Should return 2 '''
        path = self.mkbwa( 'Usage: bwa 1>&2' )
        eq_( 2, self.inst.run_bwa( [path], [], [] ) )

    def test_runbwa_validoptions( self ):
        ''' Test that bwa actually gets ran '''
        path = self.mkbwa( '$1 $2 $3' )
        eq_( 0, self.inst.run_bwa( [path], ['a','b'], ['c'] ) )

    def test_run_invalid( self ):
        ''' Make sure run is working '''
        path = self.mkbwa( 'Usage: bwa 1>&2' )
        eq_( 2, BWASubTest( bwa_path=path, command='mem' ).run() )

    def test_run_outputfile( self ):
        ''' Test that output file can be specified and is used '''
        output_file = 'output.sai'
        output_text = "test output"
        path = self.mkbwa( output_text )
        bwa = BWASubTest( bwa_path=path, command='' )
        # Should write output_text to output_file
        bwa.run( output_file )
        eq_( output_text, open( output_file ).read().rstrip('\n') )

    def test_run( self ):
        ''' Make sure run is working '''
        path = self.mkbwa( 'BWA' )
        eq_( 0, BWASubTest( bwa_path=path, command='mem' ).run() )
    
class TestBWAMem( BaseBWA ):
    @classmethod
    def setUpClass( self ):
        super( TestBWAMem, self ).setUpClass()
        # Create fake indexed fa file
        self.fa = os.path.join( self.tempdir, 'input.fa' )
        self.fai = self.fa + '.bwt'
        fh = open( self.fa, 'w' )
        fh.write( '>seq1\nATGC' )
        fh.close()
        fh = open( self.fai, 'w' )
        fh.close()
        # Create non indexed fasta
        self.fa2 = os.path.join( self.tempdir, 'input2.fa' )
        shutil.copy( self.fa, self.fa2 )

    @raises( ValueError )
    def test_requiredargs_gt2( self ):
        ''' Greater than 3 args is not correct '''
        mem = BWAMem( self.fa, self.fa2, self.fa2, 'arg4', bwa_path=self.bwa_path )
        
    @raises( ValueError )
    def test_requiredargs_lt2( self ):
        ''' Needs 2 args '''
        mem = BWAMem( self.fa, bwa_path=self.bwa_path )

    @raises( ValueError )
    def test_requiredargs_firstarg_invalid( self ):
        ''' Test invalid path for database arg '''
        BWAMem( '/invalid/path/in.fa', self.fa2, bwa_path=self.bwa_path )

    @raises( ValueError )
    def test_requiredargs_firstarg_nonindexed( self ):
        ''' Test non indexed file for index arg '''
        # Should just check to make sure self.fa2 has same filename but with .bwt
        BWAMem( self.fa2, self.fa, bwa_path=self.bwa_path )

    @raises( ValueError )
    def test_requiredargs_second_invalid( self ):
        ''' Test reads file is valid path to file '''
        BWAMem( self.fa, '/invalid/path/in.fa', bwa_path=self.bwa_path )

    @raises( ValueError )
    def test_requiredargs_third_invalid( self ):
        ''' Test reads file is valid path to file '''
        BWAMem( self.fa, self.fa2, '/invalid/path/in.fa', bwa_path=self.bwa_path )

    def test_optionspassed( self ):
        ''' Make sure all options passed make it to bwa '''
        bwa = self.mkbwa( '$@' )
        mem = BWAMem( self.fa, self.fa2, t=5, a=3, bwa_path=bwa )
        mem.run( 'output' )
        with open( 'output' ) as fh:
            result_output = fh.read().strip()
        assert '-t 5' in result_output
        assert '-a 3' in result_output
        assert self.fa + ' ' + self.fa2 in result_output

    def test_optionalthird( self ):
        ''' Test that given a correct third argument it still runs '''
        infa = ungzip( INPUT_PATH )
        bwa.index_ref( REF_PATH )
        # I don't have a mates file so just use infa again
        mem = BWAMem( REF_PATH, infa, infa, bwa_path=BWA_PATH )
        eq_( 0, mem.run() )

    def test_bwamem_run( self ):
        ''' Make sure it actually runs bwa with correct input '''
        infa = ungzip( INPUT_PATH )
        bwa.index_ref( REF_PATH )
        mem = BWAMem( REF_PATH, infa, bwa_path=BWA_PATH )
        eq_( 0, mem.run() )

    @raises(ValueError)
    def test_run_nonfastainput( self ):
        ''' Invalid fasta input file(file exists but not fasta/fastq '''
        filename = 'test.fa'
        with open( filename, 'w' ) as fh:
            fh.write( 'not a fasta' )
        bwa = BWAMem( REF_PATH, filename, bwa_path=BWA_PATH, command='mem' )

    def test_bwareturncode_count( self ):
        ''' Fixed tests for pattern matching '''
        # First item should be the sum of read counts
        tests = [
            (1, '[M::main_mem] read 1 sequences (111350 bp)...'),
            (1000, '[M::main_mem] read 100 sequences (111350 bp)...\n[M::main_mem] read 900 sequences (111350 bp)...'),
        ]
        for reads, testline in tests:
            filename = 'fasta{}.fa'.format(reads)
            util.create_fakefasta( filename, reads )
            print "Reads In File: {}".format(reads)
            print "Read Lines: {}".format(testline)
            bwa = BWAMem( self.fa, filename, bwa_path=BWA_PATH )
            eq_( 0, bwa.bwa_return_code( testline ) )

    def test_bwareturncode_nocount( self ):
        ''' Make sure if no read sequences lines exist error is returned '''
        endline = '[main] Version: {}.{}.{}-r{}'
        bwa = BWAMem( self.fa, self.fa2, bwa_path=BWA_PATH )
        eq_( 1, bwa.bwa_return_code( endline.format( 0, 7, 4, 385 ) ) )

    def test_bwareturncode_noendline( self ):
        bwa = BWAMem( self.fa, self.fa2, bwa_path=BWA_PATH )
        eq_( 1, bwa.bwa_return_code( '[main] Version: 0.7.4-r385' ) )

class TestBWAIndex( BaseBWA ):
    @raises( ValueError )
    def test_nonexistfasta( self ):
        path = '/does/not/exist.fa'
        BWAIndex( path, bwa_path=BWA_PATH ).run()

    @raises( ValueError )
    def test_existfilenotfasta( self ):
        path = 'ref.fa'
        with open( path, 'w' ) as fh:
            fh.write( 'not a fasta' )
        BWAIndex( path, bwa_path=BWA_PATH ).run()

    @raises( ValueError )
    def test_nofastaarggiven( self ):
        BWAIndex( bwa_path=BWA_PATH ).run()

    def test_validfasta( self ):
        ref = 'ref.fa'
        # Make copy in tempdir
        shutil.copy( REF_PATH, 'ref.fa' )
        bwa = BWAIndex( ref, bwa_path=BWA_PATH )
        eq_( 0, bwa.run() )
        indexes = glob.glob( ref + '.*' ) 
        assert indexes != []

    def test_nooutputfile( self ):
        ''' Ensure no output file is created '''
        os.mkdir( 'dir1' )
        os.chdir( 'dir1' )
        ref = 'ref.fa'
        # Make copy in tempdir
        shutil.copy( REF_PATH, 'ref.fa' )
        bwa = BWAIndex( ref, bwa_path=BWA_PATH )
        eq_( 0, bwa.run() )
        assert not os.path.exists( 'bwa.sai' )
        os.chdir( '..' )

class TestCompileReads( BaseBWA ):
    @classmethod
    def setUpClass( self ):
        super( TestCompileReads, self ).setUpClass()
        self.fastq, self.sff, self.ref = util.unpack_files()

    def tearDown( self ):
        for fq in glob.glob( '*.fastq' ):
            os.unlink( fq )

    def _isfastq( self, filepath ):
        ''' Simple verification that file is a fastq '''
        with open( filepath ) as fh:
            first = fh.read(1)
            eq_( '@', first )

    def test_noreads( self ):
        ''' Directory empty '''
        os.mkdir( 'emptydir' )
        outfile = bwa.compile_reads( 'emptydir' )
        eq_( [], outfile )

    def test_outputbasename( self ):
        ''' Test that the outputfile parameter works with just a basename '''
        os.mkdir( 'basename' )
        os.symlink( self.sff, os.path.join( 'basename', 'sff1.sff' ) )
        outfile = bwa.compile_reads( 'basename', 'outputfile.fastq' )
        eq_( 'outputfile.fastq', outfile )
        eq_( seqio.reads_in_file( self.sff ), seqio.reads_in_file( outfile ) )

    def test_outputabsrelpath( self ):
        ''' Test that the outputfile parameter works with a relative or abs path '''
        os.mkdir( 'absrel' )
        os.symlink( self.sff, os.path.join( 'absrel', 'sff1.sff' ) )
        outfile = bwa.compile_reads( 'absrel', 'absrel/outputfile.fastq' )
        eq_( 'absrel/outputfile.fastq', outfile )
        eq_( seqio.reads_in_file( self.sff ), seqio.reads_in_file( 'absrel/outputfile.fastq' ) )

    def test_paramsinglesff( self ):
        ''' Make sure single sff works '''
        outfile = bwa.compile_reads( self.sff )
        self._isfastq( outfile )

    def test_paramsinglefastq( self ):
        ''' Input fastq should not be changed '''
        stat = os.stat( self.fastq )
        outfile = bwa.compile_reads( self.fastq )
        eq_( stat.st_mtime, os.stat( self.fastq ).st_mtime )
        eq_( stat.st_ino, os.stat( outfile ).st_ino )

    def test_paramdirsffonly( self ):
        ''' Directory of only sff files '''
        os.mkdir( 'sffs' )
        os.symlink( self.sff, os.path.join( 'sffs', 'sff1.sff' ) )
        os.symlink( self.sff, os.path.join( 'sffs', 'sff2.sff' ) )
        outfile = bwa.compile_reads( 'sffs' )
        expected_readcount = seqio.reads_in_file( self.sff ) * 2

        self._isfastq( outfile )
        eq_( expected_readcount, seqio.reads_in_file( outfile ) )

    def test_paramdirfastqonly( self ):
        ''' Directory of only fastq '''
        os.mkdir( 'fastq' )
        os.symlink( self.fastq, os.path.join( 'fastq', 'fq1.fastq' ) )
        os.symlink( self.fastq, os.path.join( 'fastq', 'fq2.fastq' ) )
        outfile = bwa.compile_reads( 'fastq' )
        expected_readcount = seqio.reads_in_file( self.sff ) * 2

        self._isfastq( outfile )
        eq_( expected_readcount, seqio.reads_in_file( outfile ) )

    def test_paramdirmixed( self ):
        ''' Directory of sff & fastq '''
        os.mkdir( 'mixed' )
        os.symlink( self.sff, os.path.join( 'mixed', 'sff1.sff' ) )
        os.symlink( self.fastq, os.path.join( 'mixed', 'fq1.fastq' ) )
        outfile = bwa.compile_reads( 'mixed' )
        expected_readcount = seqio.reads_in_file( self.sff ) + \
            seqio.reads_in_file( self.fastq )

        self._isfastq( outfile )
        eq_( expected_readcount, seqio.reads_in_file( outfile ) )

    def test_targetbug1_1( self ):
        '''
            Targets a bug where compile_reads would generate 2 identical
            fastq files(reads.fastq and sff.fastq)
        '''
        os.mkdir( 'sffs1' )
        os.symlink( self.sff, os.path.join( 'sffs1', 'sff1.sff' ) )
        os.symlink( self.sff, os.path.join( 'sffs1', 'sff2.sff' ) )
        outfile = bwa.compile_reads( 'sffs1' )
        fastqs = glob.glob( '*.fastq' )
        assert len( fastqs ) == 1, fastqs

    def test_targetbug1_2( self ):
        '''
            Targets a bug where compile_reads would generate 2 identical
            fastq files(reads.fastq and sff.fastq)
        '''
        os.mkdir( 'sffs2' )
        os.symlink( self.sff, os.path.join( 'sffs2', 'sff1.sff' ) )
        os.symlink( self.sff, os.path.join( 'sffs2', 'fq1.fastq' ) )
        outfile = bwa.compile_reads( 'sffs2' )
        fastqs = glob.glob( '*.fastq' )
        assert len( fastqs ) == 1, fastqs

    def test_targetbug1_3( self ):
        '''
            Targets a bug where compile_reads would generate 2 identical
            fastq files(reads.fastq and sff.fastq)
        '''
        os.mkdir( 'sffs3' )
        sffpth = os.path.join( 'sffs3', 'sff1.sff' )
        os.symlink( self.sff, sffpth )
        outfile = bwa.compile_reads( sffpth )
        fastqs = glob.glob( '*.fastq' )
        assert len( fastqs ) == 1, fastqs
