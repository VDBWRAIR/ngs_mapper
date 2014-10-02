import subprocess
from os.path import dirname, basename, join, abspath, isdir, isabs
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
    # Ensure dstprefix is absolute path
    if not isabs(dstprefix):
        dstprefix = abspath(dstprefix)
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
        # Relative path to compiled bwa
        bwa_path = abspath('bwa')
        # popd
        os.chdir(curdir)
        # Path bwa will be copied to
        dst_bwa_path = join(dstprefix,'bin','bwa')
        # Copy bwa executable to dstprefix/bin
        if not isdir(dirname(dst_bwa_path)):
            os.makedirs(dirname(dst_bwa_path))
        shutil.copy2(bwa_path,dst_bwa_path)
