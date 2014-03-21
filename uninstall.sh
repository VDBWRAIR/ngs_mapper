#!/usr/bin/env bash

# Where is this script?
THIS=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)

# Get virtpath
. ${THIS}/.virtpath
rm ${THIS}/.virtpath

# Remove virtual env
echo "Removing ${virtpath}"
rm -rf ${virtpath}
echo "Removing build and dist python package directories from ${THIS}/dependencies"
find ${THIS}/dependencies -type d -name build -o -name dist -exec rm -rf {} \; 2>/dev/null
# Remove this file too
rm dependencies/matplotlib-1.3.1/lib/matplotlib/mpl-data/matplotlibrc

echo "miseqpipeline is now uninstalled"
if [ ! -z "${VIRTUAL_ENV}" ]
then
    echo "!! Make sure you run deactivate to clear the virtual env from your shell !!"
fi

