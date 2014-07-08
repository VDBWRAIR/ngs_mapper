#!/bin/bash

# This gives us the current directory that this script is in
THIS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Grab the current directory the user is in
CWD=$(pwd)
# Ensure virtpath is set
virtpath=$1
if [ -z "${virtpath}" ]
then
    # Virtualenvironment path
    echo "No virtpath given so using ${THIS}/.miseqpipeline"
    virtpath=${THIS}/.miseqpipeline
fi
# For the uninstaller
echo "virtpath=${virtpath}" > ${THIS}/.virtpath
# Where to store our goodies
binpath=${virtpath}/bin
# Where are all the source files
deppath=${THIS}/dependencies
# Where to store manpages
manpath=${virtpath}/man1

function pyinstall() {
    dirpath=$1
    oldpath=$(pwd)
    cd $dirpath
    rm -rf build dist
    python setup.py install
    _RET=$?
    cd $oldpath
}

# Ensure correct python version
python -V 2>&1 | grep -qE '2.7.[3456789]'
if [ $? -ne 0 ]
then
    echo "Please ensure you have python 2.7.3 or greater but less than 3"
    exit 1
fi

# Check to make sure required commands are available
if [ -z "$(which convert)" ]
then
    echo "Please ensure that you have the ImageMagick packages installed and then rerun this installer."
    echo "Red Hat: yum install -y ImageMagick-c++ ImageMagick"
    echo "Ubuntu: apt-get install -y imagemagick"
    exit 1
fi

INCLUDE=$(echo | cpp -x c++ -Wp,-v 2>&1 | grep -v 'ignoring' | grep -v '^#' | grep -v '^End' | xargs)
for dev in zlib.h png.h curses.h freetype.h
do
    # Ensure dev is in include path
    if [ -z "$(find ${INCLUDE} -type f -name ${dev})" ]
    then
        echo "Please ensure that the package that provides $dev is installed and then rerun this installer."
        exit 1
    fi
done

# Create the virtual environment where everything will install to
# Don't use setuptools as we will install that later
virtualenv --prompt='(miseqpipeline) ' ${virtpath}
# Activate
. ${virtpath}/bin/activate

# Make sure we are in the repository directory
cd ${THIS}

# Compile samtools if the samtools binary doesn't exist
if [ ! -e ${binpath}/samtools ]
then
    #cd ${THIS}/htslib
    #make > htslib.make.log 2>&1
    cd ${deppath}/samtools
    make > samtools.make.log 2>&1
    if [ $? -ne 0 ]
    then
        echo "Samtools failed to compile. Please check the ${deppath}/samtools.make.log for more details."
        exit 1
    fi
    ln -s $(pwd)/samtools ${binpath}/samtools
fi

# Compile bwa if the bwa binary doesn't exist
if [ ! -e ${binpath}/bwa ]
then
    cd ${deppath}/bwa
    make > bwa.make.log 2>&1
    if [ $? -ne 0 ]
    then
        echo "bwa failed to compile. Please check the ${deppath}/bwa.make.log for more details."
        exit 1
    fi
    ln -s $(pwd)/bwa ${binpath}/bwa
fi

# Some manpage setup
# First cleanse the manpath dir
rm -rf ${manpath}
mkdir ${manpath}
# Find all the actual manpages and link them into the man1 directory
find . -type f -name '*.1' | while read f
do
    # Manpages start with .TH
    head -1 "$f" | grep -q '^.TH'
    if [ $? -eq 0 ]
    then
        path_to="$(cd $(dirname "$f") && pwd)/$(basename "$f")"
        ln -s "$path_to" "${manpath}/$(basename "$f")"
    fi
done

# Install all python packages(the ordering is important because of dependencies on each other)
package_list=( PyVCF numpy nose pyparsing tornado six python-dateutil matplotlib biopython pyBWA mock cutadapt )
cd ${deppath}
for package in ${package_list[@]}
do
    pdir=$(echo ${package}*)
    echo "Installing ${pdir}"
    pyinstall ${pdir} > ${pdir}/${package}.install.log 2>&1
    if [ $_RET -ne 0 ]
    then
        echo "${package} failed to install. Please check ${THIS}/dependencies/${pdir}/${package}.install.log for more details."
        exit 1
    fi
done

# Install Trimmomatic into $virtpath/lib
TRIMOPATH=${deppath}/Trimmo*
cp -R ${TRIMOPATH} ${virtpath}/lib/

# Symlink all of our goodies into venv bin
find ${THIS} -maxdepth 1 -type f -perm /u=x,g=x,o=x | while read f
do
    echo "Installing $f to $binpath"
    ln -s ${f} ${binpath}/
done

echo "miseqpipeline is now installed"
