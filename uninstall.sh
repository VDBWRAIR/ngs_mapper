#!/usr/bin/env bash

# Where is this script?
THIS=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)

# Get virtpath from setup
virtpath=$(awk -F'=' '/^virtpath/ {print $2}' setup)

# Remove virtual env
rm -rf ${virtpath}
find ${THIS}/dependencies -type d -name build -o -name dist -exec rm -rf {} \; 2>/dev/null

echo "miseqpipeline is now uninstalled"
if [ ! -z "${VIRTUAL_ENV}" ]
then
    echo "!! Make sure you run deactivate to clear the virtual env from your shell !!"
fi

