#!/bin/bash

# We actually want the full path to scripts
scripts=$(cd $(dirname $0) && pwd)

# Where are the read files at?
reads_by_sample=$1
if [ -z "$reads_by_sample" ]
then
    echo "Please give me the location to the directory that holds all reads by sample"
    exit 1
fi

# Where is the sample/reference mapping file that the user gave us
sample_ref_map_file=$2
if [ -z "$sample_ref_map_file" ]
then
    echo "Give me a file that is space delimited with samplename reference one per line"
    exit 1
fi

# How many CPUS does the computer have(Some ugly stuff here but it works)
CPUS=$(grep '^processor' /proc/cpuinfo | wc -l)
# If CPUS comes back empty just use 1
if [ -z "$CPUS" ]
then
    CPUS=1
fi
echo "Using $CPUS cpus"

# Where to place all analysis output directories
PROJDIR=Projects

# Loop over each of our samplenames and references
# File needs to be either space or tab delimeted
# Samplename ReferencePath
# Excluding any lines that begin with a #
mkdir -p $PROJDIR

# Run in parallel
# Ignore lines beginning with #
# Remove windows newlines
grep -v '^#' $sample_ref_map_file | sed 's/\r//' | while read sample reference
do
    # Make sure that sample was set
    if [ -z "${sample}" ]
    then
        echo "${sample_ref_map_file} must have an incorrect line. Could not read samplename" >&2
        continue
    fi
    # Make sure that reference was set
    if [ -z "${reference}" ]
    then
        echo "${sample_ref_map_file} must have an incorrect line. Could not read reference" >&2
        continue
    elif [ ! -f "${reference}" ]
    then
        echo "${reference} is not a file that can be read. Skipping ${sample}" >&2
        continue
    fi
    # If the sample already has been run
    #  then don't rerun it
    if [ -d ${PROJDIR}/${sample} ]
    then
        echo "Skipping ${sample} because it already exists" >&2
        continue
    fi

    echo runsample ${reads_by_sample}/${sample} ${reference} ${sample} -od ${PROJDIR}/${sample} $RUNSAMPLEOPTIONS
done | xargs -n 6 -P $CPUS -I CMD bash -c CMD

# Graph all samples
graphs.sh -norecreate
# Softlink in consensus sequences into single folder
consensuses.sh
