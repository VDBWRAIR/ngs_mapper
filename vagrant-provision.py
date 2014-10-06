#!/usr/bin/env python

import platform
import subprocess
import os
import sys
from os.path import *
import shutil
from glob import glob
import tempfile

# Provisions pipeline into Ubuntu, CentOS or RedHat VM
# Will essentially do everything in the README.md for installation
# Then runs nosetests -v miseqpipeline at the end

# These are all the apt-get packages for ubuntu
UBUNTU_SYSTEM_PACKAGES = [
        'build-essential', 'libncurses5', 'libncurses5-dev',
        'zlib1g', 'zlib1g-dev', 'libpango1.0-0', 'libpango1.0-dev',
        'libreadline6', 'libreadline6-dev', 'openssl', 'libssl-dev',
        'unzip', 'imagemagick', 'libpng12-dev', 'default-jre',
        'git',
]

# These are all the yum packages
REDHAT_SYSTEM_PACKAGES = [
    'wget', 'ncurses', 'ncurses-devel', 'zlib',
    'zlib-devel', 'freetype', 'freetype-devel',
    'readline', 'readline-devel', 'openssl',
    'openssl-devel', 'libpng', 'libpng-devel',
    'ImageMagick', 'java-1.7.0-openjdk', 'git'
]

class NotSuperUserError(Exception): pass

def get_distribution():
    return platform.linux_distribution()

def shell_cmd( cmdstr, requireroot=False ):
    # Ensure root
    if requireroot and os.getuid() != 0:
        raise NotSuperUserError("You need to be superuser to run {0}".format(
            cmdstr
        ))
    # Run command in shell
    return subprocess.check_call(
        cmdstr, shell=True
    )

def install_redhat_packages( packages ):
    shell_cmd(
        'yum groupinstall "Development tools"',
        True
    )
    pkglist = ' '.join( packages )
    cmd = 'yum install -y ' + pkglist
    shell_cmd( cmd, True )

def install_ubuntu_packages( packages ):
    pkglist = ' '.join(packages)
    cmd = 'apt-get install -y ' + pkglist
    shell_cmd( cmd, True )

def clone_pipeline( source, dst ):
    dst = expanduser(dst)
    source = expanduser(source)
    if not isdir(dst):
        cmd = 'git clone {0} {1}'.format(source,dst)
    else:
        cmd = 'd=$(pwd); cd {0} && git stash && git pull && cd $d'.format(
            dst
        )
    shell_cmd( cmd ) 
    os.chdir(dst)

def get_python_version():
    return sys.version_info[0:3]

def install_python( version='2.7.8', installprefix='$HOME' ):
    '''
    Install python and leave python tempdir behind
    '''
    # to cd back later
    cwd = os.getcwd()
    
    # where to install
    prefix = expandvars(installprefix)
    
    # Where python will be located
    pythonexe = join(prefix,'bin','python')

    # Check to see if it is installed already
    if isfile(pythonexe):
        return pythonexe

    # Download and unpack
    shell_cmd(
        'wget --no-check-certificate https://www.python.org/ftp/python/{0}/' \
        'Python-{0}.tgz -O- | tar xzf -'.format(version)
    )

    os.chdir('Python-{0}'.format(version))
    shell_cmd(
        './configure --prefix {0}'.format(installprefix)
    )
    shell_cmd(
        'make'
    )
    shell_cmd(
        'make install'
    )

    #popd
    os.chdir(cwd)

    # Path to python executable
    return pythonexe

def create_virtualenv( venvpath='$HOME/.miseqpipeline', pythonprefix='$HOME' ):
    '''
    Unpack and install a virtualenv to venvpath
    '''
    # to cd back later
    cwd = os.getcwd()

    # Where to install path
    venvpath = expandvars(venvpath)

    # Download and unpack
    shell_cmd(
        'wget --no-check-certificate https://pypi.python.org/packages/source/v/' \
        'virtualenv/virtualenv-1.11.6.tar.gz -O- | tar xzf -'
    )

    # Create virtualenv
    pythonpath = expandvars('{0}/bin/python'.format(pythonprefix))
    shell_cmd(
        '{0} virtualenv-1.11.6/virtualenv.py' \
        ' --prompt="(miseqpipeline) " {1}'.format(pythonpath,venvpath)
    )

    # Activate virtualenv
    shell_cmd(
        '. {0}/bin/activate'.format(venvpath)
    )

    return venvpath

def install_system_packages():
    dist, vers, name = get_distribution()
    if dist == 'Ubuntu':
        install_ubuntu_packages( UBUNTU_SYSTEM_PACKAGES )
    elif dist in ('CentOS','Red Hat Enterprise Linux Workstation'):
        install_redhat_packages( REDHAT_SYSTEM_PACKAGES )
    else:
        raise ValueError("Unsupported distribution {0}".format(dist))

def run_setup( venvpath ):
    activatepath = join(venvpath,'bin','activate')
    shell_cmd(
        '. {0}; python setup.py install'.format(activatepath)
    )

def install_pipeline():
    clone_pipeline('/vagrant', '~/miseqpipeline')
    install_python()
    venvpath = create_virtualenv()
    run_setup( venvpath )

def parse_args(args=sys.argv[1:]):
    import argparse

    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group(required=True)
    
    group.add_argument(
        '--install-system-packages',
        dest='install_system',
        default=False,
        action='store_true',
        help='Just installs system packages only'
    )

    group.add_argument(
        '--install-pipeline',
        dest='install_pipeline',
        default=True,
        action='store_true',
        help='Installs pipeline'
    )

    return parser.parse_args(args)

def main( args ):
    if args.install_system:
        install_system_packages()
    else:
        install_pipeline()

if __name__ == '__main__':
    main(parse_args())
