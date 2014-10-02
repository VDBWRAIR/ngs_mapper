# Install setuptools automagically from the interwebz
from ez_setup import use_setuptools
use_setuptools()

# Use setuptools
from setuptools import setup, find_packages
from version import __version__

from glob import glob

def pip_install( pkg ):
    ''' Just run pip install pkg '''
    from subprocess import check_call, PIPE
    check_call( ['pip', 'install', pkg] )

# Numpy doesn't seem to install correctly through the install_requires section
# https://github.com/numpy/numpy/issues/2434
pip_install( 'numpy==1.8.0' )

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
    ],
    tests_require = [
    ],
    author = 'Tyghe Vallard',
    author_email = 'vallardt@gmail.com',
    description = 'Pipeline that combines sff and fastq files from multiple platforms',
    license = '',
    keywords = 'miseq iontorrent roche 454 fastq vcf',
    url = 'https://github.com/VDBWRAIR/miseqpipeline'
)
