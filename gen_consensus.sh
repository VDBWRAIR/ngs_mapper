#!/usr/bin/env bash

# Activate the precious
. /home/EIDRUdata/programs/vdbapps/bin/activate

reference=$1
if [ -z "$reference" ]
then
    echo "You did not supply me with a reference" 1>&2
    exit 1
fi

bamfile=$2
if [ -z "$bamfile" ]
then
    echo "You did not supply me with a bam file" 1>&2
    exit 1
fi

samtools mpileup -uf $reference $bamfile | bcftools view -cg - | vcfutils.pl vcf2fq
