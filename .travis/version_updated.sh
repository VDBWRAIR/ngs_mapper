#!/bin/bash

# Make sure version is changed
git diff $(basename $TRAVIS_REPO_SLUG)/__init__.py | grep -q '^.__version__' || (echo "Version was not updated"; exit 1)
