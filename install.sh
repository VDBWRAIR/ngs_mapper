#!/usr/bin/env bash

# Allow user to specify install path
# or default to current directory
INSTALL_PATH="$PWD/miniconda"
if [ ! -z "$1" ]
then
    INSTALL_PATH="$1"
fi

# Fail if any command fails
set -e
set -v

# We do this conditionally because it saves us some downloading if the
# version is the same.
if [ -z "$TRAVIS_PYTHON_VERSION" ]
then
    PYTHON_VERSION="2.7"
else
    PYTHON_VERSION="$TRAVIS_PYTHON_VERSION"
fi
if [[ "$PYTHON_VERSION" == "2.7" ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

# Install miniconda here
bash miniconda.sh -b -p "$INSTALL_PATH"

# Set path
export PATH="${INSTALL_PATH}/bin:$PATH"

# Ensure bash path search updated
hash -r

# Always say yes and don't set ps1
conda config --set always_yes yes --set changeps1 no

# Make sure conda is updated
conda update -q conda

# Useful for debugging any issues with conda
conda info -a

# Add bioconda channels(r is required for bioconda)
conda config --add channels r
conda config --add channels bioconda
conda config --add channels vdbwrair

# Install dependencies
## Conda deps first
conda install --file requirements-conda.txt
## Pip specific deps next
pip install -r requirements-pip.txt

# Install package
python setup.py install

# Tell user how to setup PATH
echo "Make sure to setup your path to include $INSTALL_PATH/bin"
echo "export PATH=$INSTALL_PATH/bin:$PATH"
