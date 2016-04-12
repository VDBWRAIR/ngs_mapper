#!/bin/bash

if [[ "${TRAVIS_PULL_REQUEST}" != "false" ]]
then
    echo "This is a pull request. Skipping release process"
    exit 0
fi

if [[ "${TRAVIS_BRANCH}" != "master" ]]
then
    echo "Not master branch. Skipping release process"
    exit 0
fi

if .travis/version_updated.sh >/dev/null 2>&1
then
    echo "Version not updated. Skipping release process"
    exit 0
fi

set -ev

git config --global user.email "builds@travis-ci.org"
git config --global user.name "Travis CI"
export GIT_TAG="v$(.travis/current_version.py)"
git tag $GIT_TAG -a -m "See [Changelog](CHANGELOG.rst) for all changes"
# Requires you first make a Personal Access tokens(https://github.com/settings/tokens)
# Then add GH_TOKEN to https://travis-ci.org/user/project/settings
export PROJECT=$(basename ${TRAVIS_REPO_SLUG})
set +v
echo "Pushing $GIT_TAG"
git push --quiet https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git $GIT_TAG
