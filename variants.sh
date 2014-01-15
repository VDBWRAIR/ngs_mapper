#!/bin/bash

# Ya, pretty ugly
CPUS=$(for pid in $(awk '/physical id/ {print $4}' /proc/cpuinfo | sort | uniq); do egrep -xA 12 "processor[[:space:]]: $pid" /proc/cpuinfo; done | awk '/cpu cores/ {print $4}' | paste -sd+ | bc)

cat samplesheet* | while read samplename ref
do
    echo miseqpipeline/varcaller.py Projects/${samplename}/${samplename}.bam ${ref} -o Projects/${samplename}/variants
done | xargs -n 5 -P $CPUS -I CMD bash -c CMD
