from glob import glob
from os.path import *
import os
import sys
import re

# Exception when no platform can be found for a read
class NoPlatformFound(Exception): pass

def filter_reads_by_platform( path, platform ):
    '''
        Filters all reads in path down to the ones that are for platform

        @returns list of reads inside of path(basename)
    '''
    # All files in path
    files = [basename(f) for f in glob( join(path,'*') ) if platform_for_read(f) == platform]
    return files

def platform_for_read( filepath ):
    '''
        Returns the platform the given read is for

        @param filepath - Path to filename
        @returns name of platform that filepath is for
    '''
    MAPPING = {
        '\S+?__[0-9]__(?:TI|RL)\d+__\d{4}_\d{2}_\d{2}__\w+.sff': 'Roche454',
        '\S+?__[0-9]__IX\d{3}__\d{4}_\d{2}_\d{2}__\w+.sff': 'IonTorrent',
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

        @returns a dictionary {'Platform': [reads for platform(basename)]
    '''
    reads_for_plat = {}
    for platform in ('Roche454', 'IonTorrent','MiSeq', 'Sanger'):
        reads = filter_reads_by_platform( path, platform )
        # Don't add them if there were none
        if reads:
            reads_for_plat[platform] = reads
    print reads_for_plat
    return reads_for_plat
