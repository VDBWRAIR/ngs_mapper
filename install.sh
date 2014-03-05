#!/bin/bash

# This gives us the current directory that this script is in
THIS="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Deactivate if possible
deactivate 2>/dev/null

# Activate vdbapps for now, but soon we will remove this in favor
# of git submodules for all the necessary apps
. /home/EIDRUdata/programs/vdbapps/bin/activate

# Grab the current directory the user is in
CWD=$(pwd)

# If any command exits then trap it and print error and exit
trap 'echo "Error running $BASH_COMMAND"; rm -rf man1; exit;' ERR SIGINT SIGTERM

# Make sure we are in the repository directory
cd ${THIS}

# Ensure zlib.h is in include path(I guess we won't assume redhat here and just do rpm -qa stuff)
if [ -z "$(find $(echo | cpp -x c++ -Wp,-v 2>&1 | grep -v 'ignoring' | grep -v '^#' | grep -v '^End' | xargs) -type f -name zlib.h)" ]
then
    echo "Please ensure that the zlib development package is installed. Probably yum install zlib-devel or apt-get install zlib-devel"
    exit 1
fi

# Now ensure submodules are setup correctly(hopefully)
# This is probably not the correct thing to do, but hacks are great yay!
# First make sure all submodules are updated and init'd
git submodule update --init
# Then just run through all of them and ensure they have HEAD checked out or something
git submodule foreach git reset --hard HEAD

# Compile samtools if the samtools binary doesn't exist
if [ ! -e ${THIS}/samtools/samtools ]
then
    cd ${THIS}/htslib
    make 2>&1 > htslib.make.log
    cd ${THIS}/samtools
    make 2>&1 > samtools.make.log
fi

# Compile bwa if the bwa binary doesn't exist
if [ ! -e ${THIS}/bwa/bwa ]
then
    cd ${THIS}/bwa
    make 2>&1 > bwa.make.log
fi

# Some manpage setup
if [ ! -d ${THIS}/man1 ]
then
    mkdir -p ${THIS}/man1
    # Find all the actual manpages and link them into the man1 directory
    find . -type f -name '*.1' | while read f
    do
        # Manpages start with .TH
        head -1 "$f" | grep -q '^.TH'
        if [ $? -eq 0 ]
        then
            path_to="$(cd $(dirname "$f") && pwd)/$(basename "$f")"
            ln -s "$path_to" "${THIS}/man1/$(basename "$f")"
        fi
    done
fi

