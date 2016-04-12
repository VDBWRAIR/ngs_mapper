#!/bin/bash

if [ "${TRAVIS_BRANCH}" != "master" ]
then
    echo "Not master branch. Skipping"
    exit 0
fi

set -e

# Make sure changelog is updated
.travis/in_commit.sh CHANGELOG.rst
