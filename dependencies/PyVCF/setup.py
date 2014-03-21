from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
    CYTHON = True
except:
    CYTHON = False

requires = []

# python 2.6 does not have argparse
try:
    import argparse
except ImportError:
    requires.append('argparse')

import collections
try:
    collections.Counter
except AttributeError:
    requires.append('counter')
try:
    collections.OrderedDict
except AttributeError:
    requires.append('ordereddict')

# get the version without an import
VERSION = "Undefined"
DOC = ""
inside_doc = False
for line in open('vcf/__init__.py'):
    if "'''" in line:
        inside_doc = not inside_doc
    if inside_doc:
        DOC += line.replace("'''", "")

    if (line.startswith('VERSION')):
        exec(line.strip())

extras = {}
if CYTHON:
    extras['cmdclass'] = {'build_ext': build_ext}
    extras['ext_modules'] = [Extension("vcf.cparse", ["vcf/cparse.pyx"])]

setup(
    name='PyVCF',
    packages=['vcf', 'vcf.test'],
    scripts=['scripts/vcf_melt', 'scripts/vcf_filter.py'],
    author='James Casbon and @jdoughertyii',
    author_email='casbon@gmail.com',
    description='Variant Call Format (VCF) parser for Python',
    long_description=DOC,
    test_suite='vcf.test.test_vcf.suite',
    install_requires=['distribute'],
    requires=requires,
    entry_points = {
        'vcf.filters': [
            'site_quality = vcf.filters:SiteQuality',
            'vgq = vcf.filters:VariantGenotypeQuality',
            'eb = vcf.filters:ErrorBiasFilter',
            'dps = vcf.filters:DepthPerSample',
            'avg-dps = vcf.filters:AvgDepthPerSample',
            'snp-only = vcf.filters:SnpOnly',
        ]
    },
    url='https://github.com/jamescasbon/PyVCF',
    version=VERSION,
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
      ],
    keywords='bioinformatics',
    use_2to3=True,
    include_package_data=True,
    package_data = {
        '': ['*.vcf', '*.gz', '*.tbi'],
        },
    **extras
)
