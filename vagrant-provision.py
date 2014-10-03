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
    return subprocess.check_output(
        cmdstr, shell=True
    )

def install_redhat_packages():
    shell_cmd(
        'yum groupinstall "Development tools"',
        True
    )
    shell_cmd(
        'yum install -y wget ncurses{,-devel} zlib{,-devel} ' \
        'freetype{,-devel} readline{,-devel} openssl{,-devel} ' \
        'libpng{,-devel} ImageMagick java-1.7.0-openjdk git',
        True
    )

def install_ubuntu_packages():
    shell_cmd(
        'apt-get install -y build-essential libncurses5{,-dev} ' \
        'zlib1g{,-dev} libpango1.0-{0,dev} libreadline6{,-dev} ' \
        'openssl libssl-dev unzip imagemagick libpng12-dev default-jre git',
        True
    )

def clone_pipeline():
    shell_cmd(
        'git clone https://github.com/VDBWRAIR/miseqpipeline.git'
    )
    os.chdir('miseqpipeline')

def install_python( version='2.7.8', installprefix='$HOME' ):
    '''
    Install python and leave python tempdir behind
    '''
    # to cd back later
    cwd = os.getcwd()
    
    # where to install
    prefix = expandvars(installprefix)

    # Download and unpack
    shell_cmd(
        'wget --no-check-certificate https://www.python.org/ftp/python/{0}/' \
        'Python-{0}.tgz -O- | tar xzf -'.format(version)
    )

    os.chdir('Python-{0}'.format(version))
    shell_cmd(
        'make'
    )
    shell_cmd(
        'make install'
    )

    #popd
    os.chdir(cwd)

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

def install_system_packages():
    dist, vers, name = get_distribution()
    if dist == 'Ubuntu':
        install_ubuntu_packages()
    elif dist in ('CentOS','Red Hat Enterprise Linux Workstation'):
        install_redhat_packages()
    else:
        raise ValueError("Unsupported distribution {0}".format(dist))

def run_setup():
    shell_cmd(
        'python setup.py install'
    )

def install_pipeline():
    clone_pipeline()
    install_python()
    create_virtualenv()
    run_setup()

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
