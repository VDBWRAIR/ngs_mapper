#!/bin/bash

if [ "${TRAVIS_BRANCH}" != "develop" ]
then
    echo "Not develop branch. Skipping"
    exit 0
fi

set -e

.travis/in_commit.sh CHANGELOG.rst docs/
