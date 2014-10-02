import subprocess
from os.path import dirname, basename, join, abspath, isdir, isabs, exists
import shutil
import os

import tempdir

def install_bwa( source, gitsha, dstprefix, tdir=None ):
    '''
    Uses git source path to any git repo(local or http[s])
    Installs the bwa executable into dstprefix/bin
    Uses tempdir location as the temporary path to build with which will be removed
    after its usage
    '''
    copypaths = ['bwa']
    clone_checkout_make_copy(source, gitsha, dstprefix, tdir=tdir, copypaths=copypaths )

def verify_bwa_install( dstprefix ):
    '''
    Make sure that bwa is installed correctly into dstprefix
    Just checks that bin/bwa exists and is executable by current user
    '''
    expath = [('bin/bwa',os.X_OK)]
    return [] == prefix_has_files(dstprefix, expath)

def prefix_has_files( prefix, checkfiles ):
    '''
    Ensures that checkfiles are all in prefix
    '''
    missing = []
    # pushd
    cwd = os.getcwd()
    os.chdir(prefix)
    # Ensure each path exists
    for p in checkfiles:
        if isinstance(p,tuple):
            accessmode = p[1]
            filee = p[0]
        else:
            filee = p
            accessmode = os.F_OK
        if not os.access(filee, accessmode):
            missing.append(filee)
    # popd
    os.chdir(cwd)
    return missing

def install_samtools( source, gitsha, dstprefix, tdir=None ):
    '''
    Uses git source path to any git repo(local or remote)
    Installs samtools and bcftools into dstprefix/bin
    Uses tdir as temp directory for installation
    '''
    copypaths = [
        'samtools',
        'bcftools/bcftools'
    ]
    clone_checkout_make_copy(source, gitsha, dstprefix, tdir=tdir, copypaths=copypaths)

def verify_samtools_install( dstprefix ):
    '''
    Ensure that dstprefix/bin contains samtools and bcftools and both are executable
    '''
    exepaths = [
        ('bin/samtools', os.X_OK),
        ('bin/bcftools', os.X_OK),
    ]
    return [] == prefix_has_files(dstprefix, exepaths)

def clone_checkout_make_copy( source, gitsha, dstprefix, copypaths=[], tdir=None ):
    '''
    Enters temporary directory
    Clones source
    Enters basename(source)
    Checks out gitsha
    runs make
    Copies all copypaths(these are relative to cloned repo path) into dstprefix/bin
    Removes temporary directory
    '''
    # Ensure dstprefix/bin is absolute path
    if not isabs(dstprefix):
        dstprefixbin = join(abspath(dstprefix),'bin')
    else:
        dstprefixbin = join(dstprefix,'bin')
    # Ensure dstprefixbin exists
    if not isdir(dstprefixbin):
        os.makedirs(dstprefixbin)
    # Run in a tempdir that gets auto cleaned up after
    with tempdir.in_tempdir(basedir=tdir) as tdir:
        # Clone the source path
        cmd = ['git','clone',source]
        subprocess.call(cmd)
        # Cloned dir name
        cloneddir = basename(source)
        # pushd
        curdir = os.getcwd()
        # Enter cloned dir
        os.chdir(cloneddir)
        # Checkout version
        cmd = ['git','checkout',gitsha]
        subprocess.call(cmd)
        # Compile
        cmd = ['make']
        subprocess.call(cmd)
        for copypath in copypaths:
            # Copy bwa executable into dstprefix/bin
            shutil.copy2(copypath, dstprefixbin)
        # popd
        os.chdir(curdir)
