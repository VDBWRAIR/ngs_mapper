#!/bin/bash

# We actually want the full path to scripts
scripts=$(cd $(dirname $0) && pwd)
# Where be the miseq run data
rawdata="/home/EIDRUdata/NGSData/RawData/MiSeq"
# What MiSeq run
miseq_run="131210_M02261_0004_000000000-A6FWH"

# Where is the sample/reference mapping file that the user gave us
sample_ref_map_file=$1
if [ -z "$sample_ref_map_file" ]
then
    echo "Give me a file that is space delimited with samplename reference one per line"
    exit 1
fi

# Loop over each of our samplenames and references
# File needs to be either space or tab delimeted
# Samplename ReferencePath
# Excluding any lines that begin with a #
grep -v '^#' $sample_ref_map_file | while read sample reference
do
    # Make sure that sample was set
    if [ -z "${sample}" ]
    then
        echo "${sample_ref_map_file} must have an incorrect line. Could not read samplename"
        continue
    fi
    # Make sure that sample was set
    if [ -z "${reference}" ]
    then
        echo "${sample_ref_map_file} must have an incorrect line. Could not read reference"
        continue
    elif [ ! -f "${reference}" ]
    then
        echo "${reference} is not a file that can be read. Skipping ${sample}"
        continue
    fi
    # If the sample already has been run
    #  then don't rerun it
    if [ -d ${sample} ]
    then
        echo "Skipping ${sample} because it already exists"
        continue
    fi

    # Otherwise LETS DO THIS
    mkdir $sample
    # Get into the samplename directory
    pushd $sample > /dev/null
    # Map the samplename
    ${scripts}/run_bwa_on_samplename.sh $sample $reference | tee -a bwa.log
    ret=$?
    # Detect if bwa didn't run correctly
    if [ $ret -ne 0 ]
    then
        echo "bwa did not run successfully for ${sample}"
    else
        ${scripts}/gen_consensus.sh $reference ${sample}.bam > ${sample}.consensus.fastq
    fi
    # Get out of the sample directory
    popd > /dev/null
done
