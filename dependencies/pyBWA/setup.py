import os
import os.path
from distutils.core import setup

from fnmatch import fnmatch
import subprocess
import sys
import glob
from bwa.install import install_bwa

# Version file to store version information in
ver_file = os.path.join( 'bwa', '_version.py' )

# The major.minor version number
# Set to 0 here as we set/read it later on
__version__ = 0

# Install bwa into bin directory so it will be copied with all of the other
# scripts inside of bin
install_bwa( 'bin/' )

# Utility function to read the README file.
# Used for the long_description. It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def scripts( ):
    return [os.path.join( 'bin', f ) for f in os.listdir( 'bin' ) \
        if not fnmatch( f, '*.swp' ) and not fnmatch( f, '*.pyc' )]

def git_branch():
    ''' Return the current checked out branch name '''
    try:
        output = subprocess.check_output( ['git', 'branch'] ).splitlines()
    except:
        print "unable to get git branch"
        return ""

    # Get the line that the astriks is in
    branch = [x for x in output if '*' in x][0]
    branch = branch.replace( '*', '' ).strip()
    # Only return branches other than master
    if branch != 'master':
        return branch
    else:
        return ''

def set_version():
    ''' Sets the version using the current tag and revision in the git repo '''
    if os.path.isdir(".git"):
        try:
            p = subprocess.Popen(["git", "describe", "--tags", "--always"], stdout=subprocess.PIPE)
        except EnvironmentError:
            print "unable to run git"
            return
        stdout = p.communicate()[0]
        if p.returncode != 0:
            print "unable to run git"
            return

        # Full version string
        ver = stdout.strip()
        branch = git_branch()
        if branch:
            ver += '.' + git_branch()
    else:
        # Hack in case no .git dir'
        ver = '0.0.0'

    with open( ver_file, 'w' ) as fh:
        global __version__
        __version__ = ver
        fh.write( "__version__ = '%s'\n" % __version__ )

    return True

# Setup
if set_version() is None:
    sys.exit( -1 )

setup(
    name = "pyBWA",
    version = __version__,
    author = "Tyghe Vallard",
    author_email = "vallardt@gmail.com",
    description = ("Python wrapper for bwa mapper"),
    keywords = "bwa walter reed research python library",
    url = "https://github.com/VDBWRAIR/pyBWA",
    packages = [
        'bwa',
    ],
    scripts = scripts(),
    data_files = [
    ],
    requires = [
        'biopython'
    ],
)

