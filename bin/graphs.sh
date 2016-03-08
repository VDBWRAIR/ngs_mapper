#!/bin/bash

# Fail if anything fails
set -e

# ngs_mapper dir
THIS=$(cd $(dirname $0) && pwd)

if [ ! -d Projects ]
then
    echo "No Projects directory"
    exit 1
fi

# Ya, pretty ugly
CPUS=$(grep '^processor' /proc/cpuinfo | wc -l)
if [ -z "$CPUS" ]
then
    CPUS=1
fi
echo "Using $CPUS cpus"

for p in Projects/*
do
    if [ "$1" == "-norecreate" ]
    then
        echo time ${THIS}/graphsample ${p}/$(basename $p).bam -od $p -qualdepth ${p}/$(basename $p).bam.qualdepth.json
    else
        echo time ${THIS}/graphsample ${p}/$(basename $p).bam -od $p
    fi
done | xargs -n 5 -P $CPUS -I CMD bash -c CMD

# Graph the mapped and unmapped reads
graph_mapunmap Projects/*/*.qualdepth.json -o MapUnmapReads.png
# Create the Sample Coverage graphic
sample_coverage Projects/* --output SampleCoverage.png
# Create one graphic for all QualDepth graphics
convert -quality 25 -compress JPEG Projects/*/*.qualdepth.png QualDepth.pdf
# Get a graphic to see how long it took to run each sample
graph_times
