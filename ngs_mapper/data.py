from glob import glob
from os.path import *
import os
import sys
import re
import log
from Bio import SeqIO
import gzip

logger = log.setup_logger(__name__, log.get_config())

ROCHE_FILE = '\S+?(?:__[0-9]){0,1}__(?:TI|RL)\d+__\d{4}_\d{2}_\d{2}__\w+.(sff|fastq)'
'''
Matches roche sff or fastq files

sample__region__barcode__year_month_day__type.filetype
'''
ROCHE_ID = '[A-Z0-9]{14}'
'''
Matches Roche accessions which are just 14 Alpha Numeric uppercase characters

@AAAAAAAAAAAAAA
'''
IONTORRENT_FILE = '\S+?__[0-9]__IX\d{3}__\d{4}_\d{2}_\d{2}__\w+.(sff|fastq)'
'''
Matches IonTorrent file names(essentially same as roche)

sample__region__barcode__year__month__day__type.filetype
'''
IONTORRENT_ID = '[A-Z]{5}:[0-9]+:[0-9]+'
'''
Matches IonTorrent(and hopefully ionproton) accessions which should be
5 ALPHA-NUMERIC characters:Digits:Digits

@IIIII:0:0
'''
SANGER_FILE = '\S+?_[FR]\d+_\d{4}_\d{2}_\d{2}_\w+_\w+_\d{4}_[A-Z]\d{2}.(fastq|ab1)'
'''
Matches Sanger file names

sample_<F or R>number_year_month_day_virus_gene_runnumber_well.filetype
'''
SANGER_ID = '[_a-zA-Z0-9-]+'
'''
Matches sanger ids which are in no specific format other than can contain
uppercase, lowercase, digits, underscore and dash

@anything-really_can-be-used
'''
MISEQ_FILE = '\S+?_S\d+_L\d{3}_R\d_\d{3}_\d{4}_\d{2}_\d{2}.fastq'
'''
Matches MiSeq file names

sample_sampleno_lane_<R1 or R2>_set_year_month_day.fastq
'''
MISEQ_ID = '[a-zA-Z0-9_]+:[0-9]+:[a-zA-Z0-9-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+'
'''
Matches MiSeq identifiers

http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm

@EAS139:136:FC706VJ:2:5:1000:12850 1:Y:18:ATCACG
'''

FILENAME_MAPPING = {
    ROCHE_FILE: 'Roche454',
    IONTORRENT_FILE: 'IonTorrent',
    SANGER_FILE: 'Sanger',
    MISEQ_FILE: 'MiSeq'
}
''' Mapping of regular expression for a filename to the platform that it belongs to '''
READ_ID_MAPPING = (
    (ROCHE_ID, 'Roche454'),
    (IONTORRENT_ID, 'IonTorrent'),
    (MISEQ_ID, 'MiSeq'),
    (SANGER_ID, 'Sanger'),
)
''' Mapping of regular expression for a read identifier to the platform that it belongs to '''

class NoPlatformFound(Exception):
    '''Exception when no platform can be found for a read'''
    pass

def filter_reads_by_platform( path, platform ):
    '''
    Filters all reads in path down to the ones that are for platform

    @returns list of reads inside of path(basename)
    '''
    # All files in path
    files = []
    for f in glob( join( path, '*' ) ):
        try:
            pfr = platform_for_read( f )
        except NoPlatformFound as e:
            logger.warning(
                "{0} was skipped as the platform cannot be determined".format(f)
            )
            continue
        except IOError as e:
            logger.warning(
                "Error({0}) trying to determine read type for {1}." \
                " It will be skipped.".format(e,f)
            )
            continue
        if pfr == platform:
            files.append( f )
    return files

def platform_for_read( filepath ):
    '''
    Returns the platform the given read is for based on a series of regular expressions

    The filepath is opened and the first read in the file is extracted
    The identifier or accession for that read is used to identify the platform for
    which the file belongs too

    You can look through the regular expressions above for how they are determined
    
    This function is able to open .gz files

    :param str filepath: Path to read to determine platform for
    :return: Roche454|IonTorrent|Sanger|MiSeq
    :raises: NoPlatformFound
    '''
    fh, ext = file_handle(filepath)
    # fix sanger abi extension issue
    if ext == 'ab1':
        ext = 'abi'
    try:
        # First record in readfile
        try:
            logger.debug('Reading first read from {0}'.format(filepath))
            first_record = next(SeqIO.parse(fh, ext))
        except StopIteration:
            logger.debug('It appears that {0} is empty'.format(filepath))
            raise NoPlatformFound("No platform found for empty read file {0}".format(filepath))
        except ValueError as e:
            logger.debug('{0} is not a parseable type by Biopython\'s Seqio.parse: {1}'.format(filepath, str(e)))
            raise NoPlatformFound("No platform found for invalid read file {0}".format(filepath))
            
        # Find first platform that matches
        for p, plat in READ_ID_MAPPING:
            if re.match(p, first_record.id):
                return plat
        raise NoPlatformFound("No platform found for {0}".format(filepath))
    finally:
        fh.close()

def file_handle(filepath):
    '''
    Return a normalized opened handle for a filepath as well as the extension

    :param str filepath: path to file
    :return: (opened filelike object, fileextension without .gz)
    '''
    root, ext = splitext(filepath)
    if ext == '.gz':
        root, ext = splitext(root)
        return (gzip.open(filepath, 'rb'), ext[1:])
    else:
        return (open(filepath, 'rb'), ext[1:])

def is_sanger_readfile( filepath ):
    '''
    Inspect the top read in the file and see if the quality encoding
    max > 40
    '''
    if not filepath.endswith( '.fastq' ):
        return False

    reads = SeqIO.parse( filepath, 'fastq' )
    r1 = next( reads )._per_letter_annotations['phred_quality']

    return max(r1) > 40

def reads_by_plat( path ):
    '''
    Returns a mapping of Platform => [reads for platform] for the given path
    If there are paired end read files they will be paired inside of a tuple in the list(MiSeq only now)

    @returns a dictionary {'Platform': [reads for platform(basename)]
    '''
    reads_for_plat = {}
    for platform in ('Roche454','IonTorrent','MiSeq', 'Sanger'):
        reads = filter_reads_by_platform( path, platform )
        # Don't add them if there were none
        if reads:
            reads_for_plat[platform] = pair_reads( reads )
    return reads_for_plat

def pair_reads( readlist ):
    '''
        Pairs paired-end read files inside of readlist into 2 item tuples
        Specific to MiSeq filenames at this point

        @param readlist - List of read file paths
        
        @returns the readlist with any reads that are paired end inside of a 2 item tuple with the first(forward) in [0]
            and the second(reverse) in [1]
    '''
    logger.debug("Attempting to pair readlist {0}".format(readlist))
    paired_reads = []
    skiplist = []
    for i in range(len(readlist)):
        logger.info("paired_reads: {0}".format(paired_reads))
        logger.debug("i: {0}".format(i))
        # Don't do anything with reads that have been mated already
        if not readlist[i] in skiplist:
            # Is there a mate file?
            index = find_mate( readlist[i], readlist )
            logger.debug("Index in readlist for {0} is {1}".format(readlist[i],index))
            # No mate found so just add it to list
            if index == -1:
                logger.debug("Index was -1 so appending readfile")
                paired_reads.append(readlist[i])
            else:
                # Append mate pair
                logger.debug("Appending mated pairs ({0},{1})".format(readlist[i],readlist[index]))
                paired_reads.append( tuple(sorted([readlist[i],readlist[index]])) )
                # Append mate to be skipped so we don't re-add it again
                skiplist.append(readlist[index])
                logger.debug("Mate list {0}".format(skiplist))
        else:
            logger.debug("Ignoring {0} because it has already been added to the mate list".format(readlist[i]))
    return paired_reads

def find_mate( filepath, readlist ):
    '''
        Finds the index of the mate file for filepath given a readlist
        Just looks for _R[12]_ and looks for identical filename but with the opposite of whatever R? is found
    
        @param filepath - Path to a read file
        @param readlist - List of paths to reads

        @returns the index in readlist for the mate for filepath or -1 if none are found
    '''
    #cp = re.compile( '(?P<samplename>\S+?)_S\d+_L\d{3}_(?P<fr>R\d)_\d{3}_\d{4}_\d{2}_\d{2}.fastq' )
    cp = re.compile('_R([12])_')
    m = cp.search( basename(filepath) )
    if not m:
        logger.debug( "Cannot find _R1_ or _R2_ in {0}".format(
            filepath
        ))
        return -1
    else:
        found_r = m.group(1)
        
        # Look for opposite index
        if found_r == '1':
            mate_r = '2'
        elif found_r == '2':
            mate_r = '1'
        else:
            logger.warning('Found {0} instead of 1 or 2 in {1}'.format(found_r, filepath))
            return -1

        matefn = filepath.replace(
            '_R{0}_'.format(found_r),
            '_R{0}_'.format(mate_r)
        )
        logger.debug('{0} found in {1} and will search for {2}'.format(found_r, filepath, matefn))
        try:
            return readlist.index(matefn)
        except ValueError as e:
            logger.debug( "Could not find mate in {0} for {1}".format(
                readlist, filepath
            ))
            return -1
