import sys
import os
import os.path
import urllib
import subprocess
import tarfile
import tempfile
import shutil
from os.path import *

def which( bin ):
    '''
        Return path to bin if it is in the system path and executable or None if it is not
    
        @param bin - Binary to attempt to find in path

        @returns the joined path from system path and executable if found. None otherwise
    '''
    for p in os.environ['PATH'].split(':'):
        # Try this path
        p = join( p, bin )
        # Is this path valid and executable?
        if os.access( p, os.X_OK ):
            return p
    return None

def get_bwa( version ):
    '''
        Download bwa from sf given the version

        @param version - the part after bwa- and before .tar.bz2
        @returns location of downloaded file
    '''
    tempd = tempfile.mkdtemp()
    dlpath = os.path.join( tempd, 'bwa-{}.tar.bz2'.format(version) )
    base_url = 'http://sourceforge.net/projects/bio-bwa/files/bwa-{}.tar.bz2/download'
    filename, info = urllib.urlretrieve( base_url.format( version ), dlpath )
    return filename

def compile_bwa( path_to_bwasource ):
    '''
        Compiles bwa source given its location

        @param path_to_bwasource - Path to directory containing bwa source code
    '''
    cmd = 'make'
    ret = subprocess.call( [cmd], cwd=path_to_bwasource )
    if ret != 0:
        raise ValueError( "Failed to compile bwa in directory {}\n".format(
            path_to_bwasource) )

def unpack_dl( archive ):
    '''
        Unpacks archive
        Autodetects gz, bz2, tgz, tar.bz2 and tar.gz files

        @param archive - Location of archive
    '''
    supported_extensions = ('gz','bz2','tgz')
    # Get archive extension
    arch, arch_ext = os.path.splitext( archive )
    if arch_ext[1:] not in supported_extensions:
        raise ValueError( "Could not determine archive type for {} Supported " \
            "extensions are {}.".format( archive, supported_extensions) )

    tf = tarfile.open( archive )
    tf.extractall( os.path.dirname( archive ) )
    tf.close()
    arch_dir = os.path.dirname( archive )
    sourced = [d for d in os.listdir( arch_dir ) \
        if os.path.isdir( os.path.join( arch_dir, d ))]
    if len( sourced ) > 1:
        raise ValueError( "Error after unpacking archive. More than 1 directory" \
            " unpacked" )
    return os.path.join( arch_dir, sourced[0] )

def install_bwa( where_to_install='/usr/local/bin', version='0.7.4', upgrade=False ):
    '''
        Downloads and installs version of bwa into where_to_install
        Just copies the executables to where_to_install after compile

        @param where_to_install - Path to install bwa compiled executables to
        @param version - version string between bwa- and .tar.bz2
        @param upgrade - True will redownload and install over the top of existing
            install
    '''
    # Don't reinstall unless told to
    if not upgrade and which( 'bwa' ):
        return
    binaries = ('bwa', 'qualfa2fq.pl', 'xa2multi.pl')
    print "Downloading bwa"
    dlfile = get_bwa( version )
    print "Unpacking bwa"
    source = unpack_dl( dlfile )
    print "Compiling bwa"
    compile_bwa( source )
    print "Installing bwa"
    for binary in binaries:
        shutil.copy( os.path.join( source, binary ), 
            os.path.join( where_to_install, binary ) )
