#!/bin/bash

# miseqpipeline dir
THIS=$(cd $(dirname $0) && pwd)

if [ ! -d Projects ]
then
    echo "No Projects directory"
    exit 1
fi

# Ya, pretty ugly
CPUS=$(for pid in $(awk '/physical id/ {print $4}' /proc/cpuinfo | sort | uniq); do egrep -xA 12 "processor[[:space:]]: $pid" /proc/cpuinfo; done | awk '/cpu cores/ {print $4}' | paste -sd+ | bc)

for p in Projects/*
do
    if [ "$1" == "-norecreate" ]
    then
        echo /usr/bin/time ${THIS}/graphsample.py ${p}/$(basename $p).bam -od $p -qualdepth ${p}/$(basename $p).bam.qualdepth.json
    else
        echo /usr/bin/time ${THIS}/graphsample.py ${p}/$(basename $p).bam -od $p
    fi
done | xargs -n 5 -P $CPUS -I CMD bash -c CMD

${THIS}/graph_mapunmap.py Projects/*/*.qualdepth.json -o MapUnmapReads.png
# Now use montage to tile all the qualdepth images into one massive one
# -geometry +1+1 does something to make the images not have borders
# -tile 1x makes them tile vertically otherwise they will come out as a matrix
montage -quality 25 -compress JPEG -geometry +1+1 -tile 1x Projects/*/*.png QualDepth.png
