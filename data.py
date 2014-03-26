from glob import glob
from os.path import *
import os
import sys
import re
import logging

logger = logging.getLogger(__name__)

# Exception when no platform can be found for a read
class NoPlatformFound(Exception): pass

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
            logger.warning( "{} was skipped as the platform cannot be determined from it's name".format(f) )
            continue
        if pfr == platform:
            files.append( f )
    #files = [f for f in glob( join(path,'*') ) if platform_for_read(f) == platform]
    return files

def platform_for_read( filepath ):
    '''
        Returns the platform the given read is for

        @param filepath - Path to filename
        @returns name of platform that filepath is for
    '''
    MAPPING = {
        '\S+?__[0-9]__(?:TI|RL)\d+__\d{4}_\d{2}_\d{2}__\w+.(sff|fastq)': 'Roche454',
        '\S+?__[0-9]__IX\d{3}__\d{4}_\d{2}_\d{2}__\w+.(sff|fastq)': 'IonTorrent',
        '\S+?_[FR]\d+_\d{4}_\d{2}_\d{2}_\w+_\w+_\d{4}_[A-Z]\d{2}.fastq': 'Sanger',
        '\S+?_S\d+_L\d{3}_R\d_\d{3}_\d{4}_\d{2}_\d{2}.fastq': 'MiSeq'
    }
    for p, plat in MAPPING.items():
        if re.match( p, basename(filepath) ):
            return plat
    raise NoPlatformFound("No platform found for {}".format(filepath))

def reads_by_plat( path ):
    '''
        Returns a mapping of Platform => [reads for platform] for the given path
        If there are paired end read files they will be paired inside of a tuple in the list(MiSeq only now)

        @returns a dictionary {'Platform': [reads for platform(basename)]
    '''
    reads_for_plat = {}
    for platform in ('Roche454', 'IonTorrent','MiSeq', 'Sanger'):
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
    paired_reads = []
    skiplist = []
    for i in range(len(readlist)):
        logger.info("paired_reads: {}".format(paired_reads))
        logger.debug("i: {}".format(i))
        # Don't do anything with reads that have been mated already
        if not readlist[i] in skiplist:
            # Is there a mate file?
            index = find_mate( readlist[i], readlist )
            logger.debug("Index in readlist for {} is {}".format(readlist[i],index))
            # No mate found so just add it to list
            if index == -1:
                logger.debug("Index was -1 so appending readfile")
                paired_reads.append(readlist[i])
            else:
                # Append mate pair
                logger.debug("Appending mated pairs ({},{})".format(readlist[i],readlist[index]))
                paired_reads.append( (readlist[i],readlist[index]) )
                # Append mate to be skipped so we don't re-add it again
                skiplist.append(readlist[index])
                logger.debug("Mate list {}".format(skiplist))
        else:
            logger.debug("Ignoring {} because it has already been added to the mate list".format(readlist[i]))
    return paired_reads

def find_mate( filepath, readlist ):
    '''
        Finds the index of the mate file for filepath given a readlist
    
        @param filepath - Path to a read file
        @param readlist - List of paths to reads

        @returns the index in readlist for the mate for filepath or -1 if none are found
    '''
    cp = re.compile( '(?P<samplename>\S+?)_S\d+_L\d{3}_(?P<fr>R\d)_\d{3}_\d{4}_\d{2}_\d{2}.fastq' )
    m = cp.match( basename(filepath) )
    if not m:
        return -1
    else:
        i = m.groupdict()
        sn = i['samplename']
        fr = i['fr']
        fri = fr[1]
        if fri == '1':
            fri = '2'
        elif fri == '2':
            fri = '1'
        else:
            return -1
        matefn = filepath.replace(fr, 'R'+fri)
        try:
            return readlist.index(matefn)
        except ValueError as e:
            return -1
