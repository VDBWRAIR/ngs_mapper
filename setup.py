from glob import glob
from os.path import join
import sys

from setuptools import setup, find_packages
import setuptools

import ngs_mapper
from ngs_mapper import util

from setuptools.command.install import install as _install
from setuptools.command.develop import develop as _develop

class install(_install):
    def run(self):
        self.run_command('build_sphinx')
        _install.run(self)

class develop(_develop):
    def run(self):
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
        'install': install,
        'develop': develop
    }
)
