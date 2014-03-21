from nose.tools import eq_, raises
from Bio import SeqIO

import shutil
import os
import os.path
import glob

import util
from .. import seqio

def is_fastq( fqpath ):
    try:
        next( SeqIO.parse( fqpath, 'fastq' ) )
        return True
    except:
        print "{} is not a fastq file".format(fqpath)
        return False

def sff_eq_fastq( sff, fastq ):
    ''' Just check read counts. sff can be a list of sff '''
    if isinstance( sff, str ):
        return reads_in( sff, 'sff' ) == reads_in( fastq, 'fastq' )
    else:
        return sum( [reads_in( s, 'sff' ) for s in sff] ) == reads_in( fastq, 'fastq' )

def reads_in( path, ftype ):
    ''' get num reads in path and ftype is fastq or sff '''
    return sum( [1 for seq in SeqIO.parse( path, ftype )] )

class SeqIOBase( util.Base ):
    @classmethod
    def setUpClass( self ):
        super( SeqIOBase, self ).setUpClass()
        self.sff_input = os.path.join( self.tempdir, 'input.sff' )
        util.ungzip( util.INPUT_SFF_PATH, self.sff_input )

    def setUp( self ):
        os.mkdir( 'tests' )
        os.chdir( 'tests' )

    def tearDown( self ):
        os.chdir( self.tempdir )
        shutil.rmtree( 'tests' )

class TestReadsInFile( SeqIOBase ):
    def readsinfiletest( self, fafile, expectedlines ):
        eq_( seqio.reads_in_file( fafile ), expectedlines )

    def test_readsinfile_fasta( self ):
        ''' Make sure simple fasta is counted correctly '''
        fh = open( 'fasta.fa', 'w' )
        fh.write( '>seq1\nABCD\n>seq2ACDE\n' )
        fh.close()
        self.readsinfiletest( 'fasta.fa', 2 )

    def test_readsinfile_fastq( self ):
        ''' Make sure a simple fastq file is counted correctly '''
        fh = open( 'fasta.fastq', 'w' )
        fh.write( '@seq1\nABCD\n+\nIIII\n@seq2\nACDE\n+\nIIII\n' )
        fh.close()
        self.readsinfiletest( 'fasta.fastq', 2 )

    def test_readsinfile_fastq2( self ):
        ''' @ symbol is part of quality and can occur at beginning of lines '''
        fh = open( 'fasta.fastq', 'w' )
        fh.write( '@seq1\nABCD\n+\n@III\n@seq2\nACDE\n+\nIIII\n' )
        fh.close()
        self.readsinfiletest( 'fasta.fastq', 2 )

    @raises(ValueError)
    def test_invalidinput( self ):
        ''' Check to make sure it detects non sff,fasta,fastq '''
        fh = open( 'fasta.fastq', 'w' )
        fh.write( 'not a valid input' )
        fh.close()
        self.readsinfiletest( 'fasta.fastq', 0 )

    @raises(ValueError)
    def test_invalidinput( self ):
        ''' Check to make sure it detects non sff,fasta,fastq '''
        fh = open( 'fasta.fastq', 'w' )
        fh.write( 'not a valid input\nnope' )
        fh.close()
        self.readsinfiletest( 'fasta.fastq', 0 )


class TestSffsToFastq( SeqIOBase ):
    @raises( ValueError )
    def test_invalid_sff( self ):
        ''' One of the sff paths in list is not valid sff path '''
        sff_list = [self.sff_input, 'invalid.sff']
        seqio.sffs_to_fastq( sff_list )
    
    @raises( ValueError )
    def test_nonlist( self ):
        ''' Param is not a list '''
        seqio.sffs_to_fastq( 'i am a string' )
        seqio.sffs_to_fastq( 1 )
        seqio.sffs_to_fastq( (1,2) )

    def test_singleitem( self ):
        ''' Param list len == 1, specify output '''
        self.runit( [self.sff_input] )

    def test_outputparam( self ):
        ''' Make sure output param works '''
        self.runit( [self.sff_input], output='myoutput.fastq' )

    @raises( ValueError )
    def test_invalidoutput( self ):
        ''' Invalid output path given '''
        self.runit( [self.sff_input], '/not/valid/path.fastq' )

    def test_multiitem( self ):
        ''' Param list len > 1, specify output '''
        shutil.copy( self.sff_input, 'input2.sff' )
        sff = ['input2.sff', self.sff_input]
        self.runit( sff )

    def runit( self, sff, output=None ):
        ''' Run sffs_to_fastq '''
        if output is None:
            path = seqio.sffs_to_fastq( sff )
        else:
            path = seqio.sffs_to_fastq( sff, output )
            assert path == output
        assert path is not None
        assert is_fastq( path )
        assert sff_eq_fastq( sff, path )

    def test_zeroitem( self ):
        ''' Param list is empty list '''
        path = seqio.sffs_to_fastq( [], output='shouldnotexist.fastq' )
        eq_( None, path )
        eq_( True, not os.path.isfile( 'shouldnotexist.fastq' ) )

class TestGetReads( SeqIOBase ):
    @raises( ValueError )
    def test_invaliddirpath( self ):
        ''' dir_path param is non existant '''
        seqio.get_reads( '/not/exist' )

    @raises( ValueError )
    def test_dirpath_notdir( self ):
        ''' dir_path is not a directory '''
        open( 'somefile', 'w' ).close()
        seqio.get_reads( 'somefile' )

    def test_direxistnofiles( self ):
        ''' No fastq or sff in dir_path but has other files '''
        os.mkdir( 'somedir' )
        open( 'somedir/somefile.fasta', 'w' ).close()
        open( 'somedir/somefile.fast', 'w' ).close()
        open( 'somedir/somefile.s', 'w' ).close()
        open( 'somedir/somefile.sffo', 'w' ).close()
        eq_( [], seqio.get_reads( 'somedir' ) ) 

    def test_emptydir( self ):
        ''' Empty directory '''
        os.mkdir( 'somedir2' )
        eq_( [], seqio.get_reads( 'somedir2' ) )

    def test_returnlist( self ):
        ''' Should return list of sff and fastq files with dir_path prefix '''
        os.mkdir( 'hasreads' )
        reads = [
            'hasreads/file1.fastq',
            'hasreads/file2.fastq',
            'hasreads/file3.sff',
            'hasreads/file4.sff',
        ]
        for read in reads:
            open( read, 'w' ).close()
        eq_( sorted(reads), sorted(seqio.get_reads( 'hasreads' )) )

class TestConcatFiles( SeqIOBase ):
    def writesomefiles( self, content, numfiles=3 ):
        ''' Write 3 files all with content in them and then return list of their paths '''
        os.mkdir( 'files' )
        os.chdir( 'files' )
        # Write a few files
        for i in range( numfiles ):
            with open( str(i), 'w' ) as fh:
                fh.write( content )
        os.chdir( '../' )
        return glob.glob( 'files/*' )
        
    @raises( ValueError )
    def test_emptyfilelist( self ):
        ''' filelist is [] '''
        seqio.concat_files( [], 'outfile.cat' )

    @raises( ValueError )
    def test_invalidfilelist( self ):
        ''' fillist is not a list '''
        seqio.concat_files( self.sff_input, 'outfile.cat' )

    @raises( ValueError )
    def test_doesnotexistsome( self ):
        ''' file in filelist does not exist but others do '''
        filelist = self.writesomefiles('test') + ['idontexist']
        seqio.concat_files( filelist, 'output' )

    @raises( ValueError )
    def test_invalidoutputfile( self ):
        ''' outputfile is not a string '''
        filelist = self.writesomefiles( 'test' )
        seqio.concat_files( filelist, ['output'] )

    @raises( seqio.EmptyFileError )
    def test_emptyoutputfile( self ):
        ''' filelist containes all empty files '''
        filelist = self.writesomefiles( '' )
        seqio.concat_files( filelist, 'output' )

    def test_emptyoutputfile_deleteafter( self ):
        ''' If output file was empty and raised error, it should clean up after itself '''
        filelist = self.writesomefiles( '' )
        try:
            seqio.concat_files( filelist, 'output' )
            assert False
        except seqio.EmptyFileError:
            assert not os.path.exists( 'output' ), 'output file still exists'

    def test_doesitwork( self ):
        ''' Does it actually concat files '''
        content = 'ContentLine1\nContentLine2\n'
        filelist = self.writesomefiles( content, 3 )
        seqio.concat_files( filelist, 'output' )
        with open( 'output' ) as fh:
            eq_( content * 3, fh.read() )

    @raises(ValueError)
    def test_inputsameoutput( self ):
        '''
            Should check that filenames are not the same
            If they are, it could result in an infinite loop that
             will consume all available free space
        '''
        filelist = self.writesomefiles( 'Text\n'*1000, 1 )
        seqio.concat_files( filelist, filelist[0] )
