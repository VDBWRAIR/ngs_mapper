#!/usr/bin/env bash

# Generate flagstats.txt inside of each project
for p in projects/*; do sn=$(basename $p); samtools flagstat ${p}/${sn}.bam > ${p}/flagstats.txt; done;

# Sample of what is returned from samtools flagstats
#334718 + 0 in total (QC-passed reads + QC-failed reads)
#0 + 0 duplicates
#329364 + 0 mapped (98.40%:-nan%)
#334718 + 0 paired in sequencing
#167369 + 0 read1
#167349 + 0 read2
#328898 + 0 properly paired (98.26%:-nan%)
#329040 + 0 with itself and mate mapped
#324 + 0 singletons (0.10%:-nan%)
#0 + 0 with mate mapped to a different chr
#0 + 0 with mate mapped to a different chr (mapQ>=5)

# This loops through all of the available statistics terms returned by samtools flagstats for each project's bam
# Runs them all through awk and counts them up
for term in 'in total' 'duplicates' ' mapped (' 'paired in sequencing' 'read1' 'read2' 'properly paired' 'with itself and mate mapped' 'singletons'
do
    echo $term
    grep "$term" projects/*/flagstats.txt | awk -F':' '{print $2}' | awk '{count+=$1} END {print count}'
done
