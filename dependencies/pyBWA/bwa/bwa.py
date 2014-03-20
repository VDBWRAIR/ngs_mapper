import logging
import re
from subprocess import Popen, PIPE, check_output
import tempfile
import os
import os.path
import sys
import glob
import fnmatch
import tempfile

try:
    from Bio import SeqIO
except ImportError:
    print "Please ensure biopython is installed"
    print "Try pip install biopython"
    sys.exit(1)

import seqio

logger = logging.getLogger( __name__ )

def compile_reads( reads, outputfile='reads.fastq' ):
    '''
        Compile all given reads from directory of reads or just return reads if it is fastq
        If reads is sff file then convert to fastq

        @param reads - Directory/file of .fastq or .sff
        @param outputfile - File path of single fastq file output
        @return fastq with all reads from reads
    '''
    if os.path.isdir( reads ):
        reads = seqio.get_reads( reads )
    elif isinstance( reads, str ):
        # Single read file given
        if os.path.splitext( reads )[1] == '.sff':
            # Just convert the single reads
            return seqio.sffs_to_fastq( [reads], outputfile )
        else:
            # Already fastq so nothing to do
            #  This is a bad assumption
            return reads
    
    # Empty read list
    if not len( reads ):
        return []

    # Get only sff files to convert
    sffs = fnmatch.filter( reads, '*.sff' )
    tmpsfffastq = None
    if len( sffs ):
        tmpsfffastq = os.path.join(
            os.path.dirname( outputfile ),
            'sff.' + os.path.basename( outputfile )
        )
        logger.info( "Concatting and Converting {} to fastq".format(sffs) )
        sfffastq = [seqio.sffs_to_fastq( sffs, tmpsfffastq )]
    else:
        sfffastq = []

    fastqs = fnmatch.filter( reads, '*.fastq' )
    # Concat fastq files and sff converted fastq files into
    #  outputfile
    converts = fastqs + sfffastq
    logger.info( "Concatting {} to {}".format(
            converts, outputfile
        )
    )
    seqio.concat_files( fastqs + sfffastq, outputfile )
    if tmpsfffastq is not None:
        os.unlink( tmpsfffastq )
    return outputfile

def compile_refs( refs ):
    '''
        Compile all given refs into a single file to be indexed

        @TODO -- Write tests

        @param refs - Directory/file of fasta formatted files
        @return path to concatted indexed reference file
    '''
    ref_files = []
    ref_extensions = ('.fa', '.fasta', '.fna', '.fas')

    if os.path.isdir( refs ):
        logger.info( "Compiling and concatting refs inside of {}".format(refs) )
        files = glob.glob( os.path.join( refs, '*' ) )
        logger.debug( "All files inside of {}: {}".format( files, refs ) )
        ref_files = [f for f in files if os.path.splitext(f)[1] in ref_extensions]
        logger.debug( "Filtering files down to only files with extensions in {}".format(ref_extensions) )
        logger.debug( "Filtered files to concat: {}".format( ref_files ) )
        try:
            seqio.concat_files( ref_files, 'reference.fa' )
        except (OSError,IOError,ValueError) as e:
            logger.error( "There was an error with the references in {}".format(refs) )
            logger.error( str( e ) )
            sys.exit(1)
        return 'reference.fa'
    else:
        return refs

def is_indexed( ref ):
    '''
        Checks to see if a given reference is indexed already

        @param - Refrence file name
        @return True if ref is indexed, False if not
    '''
    ref_ext = ('.amb', '.ann', '.bwt', '.pac', '.sa')
    ref_indexes = glob.glob( ref + '.*' )

    return len( ref_indexes ) != 0 and all( [True for index in ref_indexes if os.path.splitext( index )[1] in ref_ext] )

def index_ref( ref, bwa_path=None ):
    '''
        Indexes a given reference

        @param ref - Reference file path to index
        @param bwa_path - Optional path to bwa executable
    '''
    # Don't reindex an already indexed ref
    if is_indexed( ref ):
        return True

    if bwa_path is None:
        bwa_path = which_bwa()
        logger.debug( "BWA path not specified so using default " \
            " path {}".format( bwa_path ) )

    logger.info( "Indexing {}".format(ref) )
    try:
        ret = BWAIndex( ref, bwa_path=bwa_path ).run()
    except ValueError as e:
        logger.error( e )

    if ret != 0:
        logger.error( "Error running bwa index on {}".format( ref ) )
        return False
    else:
        logger.info( "bwa index ran on {}".format(ref) )
        return True

def which_bwa( ):
    '''
        Return output of which bwa
    '''
    return check_output( ['which', 'bwa'] ).strip()

def bwa_usage():
    '''
        Returns the output of just running bwa mem from command line
    '''
    return check_output( ['bwa', 'mem'] ).strip()

class BWA( object ):
    # Options that are required
    REQUIRED_OPTIONS = ['bwa_path', 'command']
    # regex to detect usage output
    USAGE_REGEX = re.compile( 'Usage:\s*bwa' )

    def __init__( self, *args, **kwargs ):
        '''
            Base class to run mem or aln options

            BWA Options:
                kwargs represent any option that has a dash before where the option
                    with the dash is the key and the value is the value

                args represent the required options(handled in subclasses)
                the first of these should be one of the main commands(mem, aln...)
            
            Class Options:
                bwa_path as a kwarg that specifies the path to the bwa executable
        '''
        # Save args, kwargs for parsing
        self.kwargs = kwargs
        self.args = list( args )
        # Options list
        self.options = []
        # Required options values. If you zip REQUIRED_OPTIONS and required_options_values you will get a
        #  mapping of k,v pairs
        self.required_options_values = []
        # Parse and remove required_options
        self.required_options()
        # Setup options from the rest of the kwargs
        self.compile_bwa_options()
        # This needs to be implemented in subclass
        self.required_args()

    def required_args( self ):
        ''' Sets self.args '''
        raise NotImplementedError( "This class is intended to be subclassed " \
                "and not instantiated directly" )
    
    def required_options( self ):
        '''
            Parse out REQUIRED_OPTIONS from kwargs and set them in 
            self.required_options_values
        '''
        try:
            # Build up the values in order they appear in REQUIRED_OPTIONS
            for op in self.REQUIRED_OPTIONS:
                self.required_options_values.append( self.kwargs[op] )
                # No longer need this in kwargs
                del self.kwargs[op]
        except KeyError as e:
            # Detects if a parameter is missing
            raise ValueError( "{} is a required parameter".format(op) )

    def compile_bwa_options( self ):
        '''
            Convert kwargs to options list
            Assumes REQUIRED_OPTIONS are not part of this list
            @returns list of options aka [k1, v1, k2, v2...]
        '''
        # Build up self.options from kwargs
        for op, val in self.kwargs.items():
            # None type values should be ignored
            if val is None or val == '':
                continue
            # Append dash to option
            self.options.append( '-'+op )
            # Options should all be strings(just being passed to command line anyways)
            val = str(val)
            # True false values only have option
            if val.lower() not in ('true','false'):
                self.options.append( val )

    def bwa_return_code( self, output ):
        '''
            Parse stderr output to find if it executed without errors
            Since it seems that bwa does not set return codes we have to parse
            stderr output instead

            If the following regex is found then the Usage statement was printed 
             which indicates a failure of one of the options:
                ^Usage:\s+bwa
            
            Subclasses need to implement this as well and call this but they need 
            to parse the rest of the output if this returns success in order to 
            tell if the algorithm ran correctly or not
        
            @returns 0 if no usage was found, 1 if usage was found
        '''
        # Search the output
        m = self.USAGE_REGEX.search( output )

        # If there is a match return 1
        if m:
            logger.warning( "BWA Returned Usage help instead of running. " \
                "This could indicate an error." )
            return 2
        # Otherwise return 0
        return 0

    def run( self, output_file='bwa.sai' ):
        '''
            Wrapper function to make running bwa easier so you don't have to supply
            a bunch of arguments

            @param output_file - The file path to write the sai output to
            @returns output of self.run_bwa
        '''
        return self.run_bwa( self.required_options_values, self.options, 
            self.args, output_file )

    def run_bwa( self, required_options, options_list, args_list, output_file='bwa.sai' ):
        '''
            @param required_options - Should correspond to self.REQUIRED_OPTIONS
            @param options_list - Full options for bwa as a list (ex. ['mem', '-t', '2'])
            @param args_list - Required arguments that come after options
            @param output_file - Output location for stdout

            @returns 0 for success, 2 if incorrect options

            Subclass implementation should return 1 for any other failures
        '''
        if not os.path.exists( required_options[0] ):
            raise ValueError( "{} is not a valid bwa path".format( required_options[0] ) )

        # Run bwa
        with open( output_file, 'wb' ) as fh:
            cmd = required_options + options_list + args_list
            logger.info( "Running {}".format( " ".join( cmd ) ) )
            p = Popen( cmd, stdout=fh, stderr=PIPE )

            # Get the output
            stdout, stderr = p.communicate()
        logger.debug( "STDERR: {}".format(stderr) )

        # Parse the status
        return self.bwa_return_code( stderr )

    def validate_indexed_fasta( self, fastapath ):
        '''
            Make sure fastapath is a valid path and already has an index

            @param fastapath - Path to fasta file
        '''
        if not is_indexed( fastapath ):
            raise ValueError( "{} does not have an index".format(fastapath) )
        if not os.path.exists( fastapath ):
            raise ValueError( "{} does not exist".format(fastapath) )

    def validate_input( self, inputpath ):
        if os.path.exists( inputpath ):
            try:
                seqio.seqfile_type( inputpath )
            except ValueError:
                raise ValueError( "{} is not a valid input file".format(inputpath) )
        else:
            raise ValueError( "{} is not a valid input file".format(inputpath) )

class BWAIndex( BWA ):
    def __init__( self, *args, **kwargs ):
        ''' Injects index command and runs super '''
        kwargs['command'] = 'index'
        super( BWAIndex, self ).__init__( *args, **kwargs )

    def required_args( self ):
        '''
            Index only requires an input fasta file to index
            Validate that it is an actual fasta file
        '''
        if len( self.args ) != 1:
            raise ValueError( "bwa index needs only 1 parameter" )
        self.validate_input( self.args[0] )
        if seqio.reads_in_file( self.args[0] ) == 0:
            raise ValueError( "{} is not a valid file to index".format(self.args[0]) )

    def bwa_return_code( self, stderr ):
        ''' 
            Missing file:
                [bwa_index] fail to open file 'bob'. Abort!

            bwa index runs successfully pretty much no matter what
        '''
        if '[bwa_index] fail to open file' in stderr:
            return 1
        return 0

    def run( self ):
        '''
            Call super and then remove output file
            hackish
        '''
        fd, tmpf = tempfile.mkstemp()
        ret = super( BWAIndex, self ).run( tmpf )
        os.unlink( tmpf )
        return ret

class BWAMem( BWA ):
    def __init__( self, *args, **kwargs ):
        ''' Injects mem command and runs super '''
        kwargs['command'] = 'mem'
        super( BWAMem, self ).__init__( *args, **kwargs )

    def required_args( self ):
        '''
            Mem requires 2 args
                db.prefix - Indexed reference genome
                reads.fq - fastq file to map reads from

            It also has one optional argument
                mates.fq - mates fastq file

            Ensures args are valid
        '''
        # Validate args
        # First argument has to be a valid indexed fasta
        self.validate_args( self.args )

    def validate_args( self, args ):
        if len( args ) > 3:
            raise ValueError( "Too many arguments supplied to BWAMem" )
        elif len( args ) < 2:
            raise ValueError( "Too few arguments supplied to BWAMem" )
        else:
            self.validate_indexed_fasta( self.args[0] )
            self.validate_input( self.args[1] )
            if len( args ) == 3:
                self.validate_input( self.args[2] )

    def bwa_return_code( self, output ):
        '''
            Just make sure bwa output has the following regex and make sure the read \d counts
            up to how many sequence lines there are
            Should end with Version line

            Example Line:
                [M::main_mem] read 100 sequences (111350 bp)...
                [main] Version: 0.7.4-r385
        '''
        read_line_pat = '\[M::main_mem\] read (\d+) sequences \((\d+) bp\)...'
        cpat = re.compile( read_line_pat )

        total_reads = 0
        total_bp = 0
        counts = cpat.findall( output )
        for reads, bps in counts:
            total_reads += int( reads )
            total_bp += int( bps )

        # Count num of read sequences
        expected_reads = seqio.reads_in_file( self.args[1] )
        # If mates file was given count them too
        if len( self.args ) == 3:
            expected_reads += seqio.reads_in_file( self.args[2] )

        # No lines found in input file?
        if expected_reads == 0:
            return 1

        if total_reads != expected_reads:
            logger.warning( "Expecting BWA to process {} reads but processed {}".format(expected_reads, total_reads) )
            return 1

        return super( BWAMem, self ).bwa_return_code( output )
