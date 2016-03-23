from ngs_mapper.ez_setup import use_setuptools
use_setuptools()

from glob import glob
import sys
from os.path import join, expanduser
import os
import subprocess

from setuptools import setup, find_packages
import setuptools
from setuptools.command.bdist_egg import bdist_egg as _bdist_egg
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install

import ngs_mapper
from ngs_mapper import util

class InstallSystemPackagesCommand(setuptools.Command):
    '''
    Custom setup.py install keyword to initiate system-package installation
    '''
    user_options = []
    description = 'Installs all system packages via the package manager(Requires super-user)'

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        from ngs_mapper.dependency import (
            install_system_packages,
            get_distribution_package_list,
            UserNotRootError,
            make_directory_readable
        )
        # Ensure readable/writeable for everybody since it is likely installed
        # first by root
        p = subprocess.Popen('chmod -R ugo=rwX .eggs setuptools* pipeline.log', shell=True)
        p.wait()

        try:
            system_packages = get_distribution_package_list('system_packages.lst')
            install_system_packages(system_packages)
        except UserNotRootError as e:
            print "You need to be root to install system packages"

class InstallPythonCommand(setuptools.Command):
    '''
    Wrapper around installing python into HOME prefix
    '''
    description = 'Allows the user to easily install python into a prefix'
    user_options = [
        ('prefix=', None, 'Where to install python to. Default is $HOME'),
        ('version=', None, 'What version to install. Default is 2.7.8')
    ]

    def initialize_options(self):
        self.prefix = expanduser('~/')
        self.version = '2.7.8'

    def finalize_options(self):
        pass

    def run(self):
        from ngs_mapper.dependency import install_python
        install_python( self.prefix, self.version )

class PipelineInstallCommand(_install):
    '''
    Custom install command which should install everything needed
    '''
    #user_options = []
    description = 'Installs the pipeline'

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Because setuptools is so terrible at handling this
        self.pip_install( 'requirements.txt' )
        # Install dependencies outside of python
        # May require that setup_requires has been processed
        # so has to come after _bdist_egg.run
        print "Installing ext deps"
        self._install_external_dependencies()

    def pip_install( self, reqfile ):
        ''' Just run pip install pkg '''
        from subprocess import check_call, PIPE
        cmd = ['pip', 'install', '--upgrade', '-r', reqfile]
        print ' '.join(cmd)
        check_call(cmd)

    def _install_external_dependencies(self):
        # URLs for dependencies
        bwa_url = 'https://github.com/lh3/bwa'
        samtools_url = 'https://github.com/samtools/samtools'
        trimmomatic_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip'

        # Install samtools and bwa
        from ngs_mapper.dependency import (
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
        install_samtools(samtools_url, 'ccf1da91b29b75764402e20f46ec21fc293fe5f5', prefix)
        install_trimmomatic(trimmomatic_url, libdir)

class bdist_egg(_bdist_egg):
    def run(self):
        self.run_command('install_pipeline')
        self.run_command('build_sphinx')
        _bdist_egg.run(self)

class develop(_develop):
    def run(self):
        install_pipeline = self.distribution.get_command_obj('install_pipeline')
        install_pipeline.develop = True
        self.run_command('install_pipeline')
        self.run_command('build_sphinx')
        _develop.run(self)

# Run setuptools setup
setup(
    name = ngs_mapper.__projectname__,
    version = ngs_mapper.__version__,
    packages = find_packages(),
    scripts = glob('bin/*'),
    entry_points = {
        'console_scripts': [
            'is_sanger = ngs_mapper.scripts:is_sanger',
            'convert_sangers = ngs_mapper.scripts:convert_sangers',
            'sff_to_fastq = ngs_mapper.sff_to_fastq:main',
            'convert_formats = ngs_mapper.convert_formats:main',
            'ngs_filter = ngs_mapper.nfilter:main',
            'roche_sync = ngs_mapper.roche_sync:main',
            'sample_coverage = ngs_mapper.coverage:main',
            'make_example_config = ngs_mapper.config:main',
            'base_caller = ngs_mapper.base_caller:main',
            'ion_sync = ngs_mapper.ion_sync:main',
            'fqstats = ngs_mapper.fqstats:main',
            'graph_mapunmap = ngs_mapper.graph_mapunmap:main',
            'graphsample = ngs_mapper.graphsample:main',
            'graph_times = ngs_mapper.graph_times:main',
            'miseq_sync = ngs_mapper.miseq_sync:main',
            'rename_sample = ngs_mapper.rename_sample:main',
            'run_bwa_on_samplename = ngs_mapper.run_bwa:main',
            'runsample = ngs_mapper.runsample:main',
            'sanger_sync = ngs_mapper.sanger_sync:main',
            'stats_at_refpos = ngs_mapper.stats_at_refpos:main',
            'tagreads = ngs_mapper.tagreads:main',
            'trim_reads = ngs_mapper.trim_reads:main',
            'vcf_consensus = ngs_mapper.vcf_consensus:main',
            'vcf_diff = ngs_mapper.vcf_diff:main',
        ]
    },
    setup_requires = [
        'nose',
        'mock',
        'tempdir',
        'sphinx',
        'sphinx_rtd_theme',
        'logconfig',
        'sh'
    ],
    tests_require = [
        'nose',
        'mock',
    ],
    package_data = {
        'ngs_mapper': ['config.yaml','MidParse.conf'],
    },
    author = ngs_mapper.__authors__,
    author_email = ngs_mapper.__authoremails__,
    description = ngs_mapper.__description__,
    license = ngs_mapper.__license__,
    keywords = 'miseq iontorrent roche 454 fastq vcf',
    url = ngs_mapper.__url__,
    data_files = [
    ] + util.build_datafiles(join(sys.prefix,'docs/ngs_mapper'), 'doc/build/html'),
    cmdclass = {
        'install_system_packages': InstallSystemPackagesCommand,
        'install_pipeline': PipelineInstallCommand,
        'install_python': InstallPythonCommand,
        'bdist_egg': bdist_egg,
        'develop': develop,
    },
)
