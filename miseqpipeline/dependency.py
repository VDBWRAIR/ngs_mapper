import subprocess
from os.path import (
    dirname, basename,
    join, abspath,
    isdir, isabs,
    exists, splitext
)
import shutil
import os
import urllib
import zipfile
import tarfile
from glob import glob
import platform
import json

import tempdir

def install_bwa( source, gitsha, dstprefix ):
    '''
    Uses git source path to any git repo(local or http[s])
    Installs the bwa executable into dstprefix/bin
    Uses tempdir location as the temporary path to build with which will be removed
    after its usage
    '''
    copypaths = ['bwa']
    clone_checkout_make_copy(source, gitsha, dstprefix, copypaths=copypaths )

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

def install_samtools( source, gitsha, dstprefix ):
    '''
    Uses git source path to any git repo(local or remote)
    Installs samtools and bcftools into dstprefix/bin
    Uses tdir as temp directory for installation
    '''
    copypaths = [
        'samtools',
        'bcftools/bcftools'
    ]
    clone_checkout_make_copy(source, gitsha, dstprefix, copypaths=copypaths)

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

def download_unpack( source, dest ):
    '''
    Downloads or copies source and then unpacks it into dest
    '''
    # Get the basename and extension
    base_name, ext = splitext(basename(source))

    # Ensure dest is absolute path
    if not isabs(dest):
        dest = abspath(dest)
    
    # Determine how we will read the file
    if ext == '.zip':
        print "Detected zipfile"
        opener = zipfile.ZipFile
        openmode = 'r'
    elif 'tar' in ext or 'gz' in ext or 'bz' in ext:
        print "Detected tarfile(tar, gz or bz)"
        opener = tarfile.open
        openmode = 'r:*'
    else:
        raise ValueError(
            '{0} cannot be extract as it is an unknown format'.format(
                source
        ))

    # Do it all in tempdir
    with tempdir.in_tempdir() as tdir:
        print "Entered temporary directory {0}".format(tdir)

        # Download into the fifo
        print "Downloading {0}".format(source)
        dlfile, headers = urllib.urlretrieve(source)

        # Extract from fifo into dest
        fh = opener(dlfile, openmode)
        print "Extracting all into {0}".format(dest)
        fh.extractall(dest)

def install_trimmomatic( source, dst ):
    '''
    Download and unpack trimmomatic into dstprefix/lib
    '''
    print "Unpacking trimmomatic into {0}".format(dst)
    download_unpack( source, dst )

def verify_trimmomatic( dstprefix, version='0.32' ):
    '''
    Ensures that the trimmomatic-<version>.jar is located inside of 
    dstprefix/lib/Trimmomatic-<version>
    '''
    paths = [
        ('lib/Trimmomatic-{0}/trimmomatic-{0}.jar'.format(version), os.R_OK)
    ]
    missing = prefix_has_files(dstprefix, paths)
    print missing
    return [] == missing

''' Class for when an unknown linux distribution is encountered '''
class UnknownDistributionError(Exception): pass

def get_distribution_package_manager( ):
    '''
    Return the package manager that should be used for this distribution

    Ubuntu - apt-get
    Debian - apt-get
    CentOS - yum
    Red Hat - yum

    All others will raise UnknownDistributionError
    '''
    dist, version, name = platform.linux_distribution()
    if dist in ('Ubuntu', 'debian'):
        return 'apt-get'
    elif dist.startswith('Red Hat Enterprise') or dist == 'CentOS':
        return 'yum'
    else:
        raise UnknownDistributionError(
            "{0} is an unknown distribution".format(dist)
        )

''' Class for when user needs to be root and is not '''
class UserNotRootError(Exception): pass

def install_system_packages( packagelist ):
    '''
    Uses the given packagemanager and package manager options to install the given
    package list
    '''
    # Get the package manager for current distribution
    pkgmanager = get_distribution_package_manager()

    # Ensure super-user or raise exception
    if os.getuid() != 0:
        raise UserNotRootError(
            "Using the package manager {0} requires super-user".format(
                pkgmanager
            )
        )

    # Run the package manager install without prompt(-y)
    cmd = [pkgmanager, 'install', '-y'] + packagelist
    subprocess.check_call( cmd )

''' Class for when package manager entry is missing from pkglistfile '''
class MissingPackageListEntry(Exception): pass

def get_distribution_package_list( pkglistfile ):
    '''
    Return the package list from a given pkglistfile
    Essentially just a dictionary lookup after loading a json encoded file
    '''
    pkgmanager = get_distribution_package_manager()

    # Load the json file with the package lists for each package manager
    pkglists = None
    with open(pkglistfile) as fh:
        pkglist = json.load(fh)

    # Ensure pkgmanager has an entry
    if pkgmanager not in pkglist:
        raise MissingPackageListEntry(
            "{0} does not contain a package list for package " \
            "manager {1}".format(pkglistfile,pkgmanager)
        )
    # return the list
    return pkglist[pkgmanager]

def which_newer_version( version1, version2 ):
    '''
    Return the highest version of the two
    '''
    return sorted([version1,version2])[1]

def install_python( prefix, version='2.7.8' ):
    '''
    Just runs configure, make and make install after downloading
    Python
    '''
    # Don't reinstall unless we are upgrading
    pyexepth = join(prefix,'bin','python')
    if exists(pyexepth):
        # Because python 2.6 doesn't have check_output
        p = subprocess.Popen([pyexepth,'--version'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        sout,sin = p.communicate()
        versionoutput = sout.split()[1]
        highversion = which_newer_version(version, versionoutput)
        # If the higher of the two versions is already installed then do not
        # reinstall or if the versions are identical
        if highversion != version or version == versionoutput:
            print "{0} already installed and satisfies requirements".format(versionoutput)
            return

    # pushd
    cwd = os.getcwd()

    # Download and unpack Python-version
    dl_url = 'https://www.python.org/ftp/python/{0}/Python-{0}.tgz'
    dl_url = dl_url.format(version)
    download_unpack(dl_url, os.getcwd())

    # Enter Python unpacked directory
    os.chdir('Python-{0}'.format(version))

    # Run configure
    cmd = ['./configure', '--prefix', prefix]
    subprocess.check_call(cmd)

    # Run make
    subprocess.check_call(['make'])

    # Run make install to install into prefix
    subprocess.check_call(['make','install'])

    os.chdir(cwd)
