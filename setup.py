# Install setuptools automagically from the interwebz
from ez_setup import use_setuptools
use_setuptools()

from glob import glob
import sys
from os.path import join

from setuptools import setup, find_packages
import setuptools.command.install as setuptools_install

from version import __version__

class CustomInstallCommand(setuptools_install.install):
    '''
    Custom setup.py install keyword to initiate system-package installation
    '''
    user_options = setuptools_install.install.user_options + [
        ('system-packages', None, 'Install system packages'),
    ]

    def initialize_options(self):
        setuptools_install.install.initialize_options(self)
        self.system_packages = False

    def run(self):
        if self.system_packages:
            from miseqpipeline.dependency import (
                install_system_packages,
                get_distribution_package_list,
                UserNotRootError
            )
            try:
                system_packages = get_distribution_package_list('system_packages.lst')
                install_system_packages(system_packages)
            except UserNotRootError as e:
                print "You need to be root to install system packages"
        else:
            # Numpy doesn't seem to install correctly through the install_requires section
            # https://github.com/numpy/numpy/issues/2434
            print "Installing numpy"
            self.pip_install( 'numpy==1.8.0' )
            # Run normal setuptools install
            print "Installing pipeline"
            setuptools_install.install.run(self)
            # Install dependencies outside of python
            print "Installing ext deps"
            self._install_external_dependencies()

    def _install_external_dependencies(self):
        # URLs for dependencies
        bwa_url = 'https://github.com/lh3/bwa'
        samtools_url = 'https://github.com/samtools/samtools'
        trimmomatic_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip'

        # Install samtools and bwa
        from miseqpipeline.dependency import (
                install_samtools,
                install_bwa,
                install_trimmomatic
        )

        # Prefix path for installation
        prefix = sys.prefix
        bindir = join(prefix,'bin')
        libdir = join(prefix,'lib')

        # Install all dependencies outside fo pypi
        install_bwa(bwa_url, '0.7.6a', prefix)
        install_samtools(samtools_url, '96b5f2294ac005423', prefix)
        install_trimmomatic(trimmomatic_url, libdir)

    def pip_install( self, pkg ):
        ''' Just run pip install pkg '''
        from subprocess import check_call, PIPE
        check_call( ['pip', 'install', pkg] )

# Run setuptools setup
setup(
    name = "miseqpipeline",
    version = __version__,
    packages = find_packages(),
    scripts = glob('bin/*'),
    install_requires = [
        'PyVCF==0.6.6',
        'numpy==1.8.0',
        'python-dateutil==2.1',
        'matplotlib==1.3.1',
        'biopython==1.63',
        'cutadapt==1.2.1',
        'nose',
        'mock',
        'pyBWA==v0.2.2',
        'tempdir',
    ],
    dependency_links = [
        'git+https://github.com/VDBWRAIR/pyBWA#egg=pyBWA-v0.2.2',
    ],
    setup_requires = [
        'tempdir'
    ],
    tests_require = [
    ],
    author = 'Tyghe Vallard',
    author_email = 'vallardt@gmail.com',
    description = 'Pipeline that combines sff and fastq files from multiple platforms',
    license = '',
    keywords = 'miseq iontorrent roche 454 fastq vcf',
    url = 'https://github.com/VDBWRAIR/miseqpipeline',
    cmdclass = {
        'install': CustomInstallCommand,
    },
)
