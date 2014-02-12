#!/bin/bash

# Path to miseqpipeline dir
THIS=$(cd $(dirname $0) && pwd)

if [ ! -d Projects ]
then
    echo "There is no Projects directory"
    exit 1
fi

if [ ! -e samplesheet* ]
then
    echo "There is no samplesheet"
    exit 1
fi

# Ya, pretty ugly
CPUS=$(for pid in $(awk '/physical id/ {print $4}' /proc/cpuinfo | sort | uniq); do egrep -xA 12 "processor[[:space:]]: $pid" /proc/cpuinfo; done | awk '/cpu cores/ {print $4}' | paste -sd+ | bc)

cat samplesheet* | while read samplename ref
do
    echo ${THIS}/varcaller.py Projects/${samplename}/${samplename}.bam ${ref} -o Projects/${samplename}/variants
done | xargs -n 5 -P $CPUS -I CMD bash -c CMD
