#!/bin/bash

for p in Projects/*
do
    if [ "$1" == "-norecreate" ]
    then
        echo /usr/bin/time miseqpipeline/graphsample.py ${p}/$(basename $p).bam -od $p -qualdepth ${p}/$(basename $p).bam.qualdepth.json
    else
        echo /usr/bin/time miseqpipeline/graphsample.py ${p}/$(basename $p).bam -od $p
    fi
done | xargs -n 5 -P 12 -I CMD bash -c CMD

convert -quality 25 -compress JPEG Projects/*/*.png QualDepth.pdf
