#!/usr/bin/env python
# just print current value of __version__ in __init__.py
import os
from os.path import join, basename, dirname

THIS = dirname(__file__)
TRAVIS_REPO_SLUG = os.environ['TRAVIS_REPO_SLUG']
project_name = basename(TRAVIS_REPO_SLUG)

execfile(join(THIS, '../', project_name, '__init__.py'))
print __version__
