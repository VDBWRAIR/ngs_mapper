#!/bin/bash

# Remove old consensus directory
if [ -d vcf_consensus ]
then
    rm -rf vcf_consensus
fi
# Make consensus dir
mkdir -p vcf_consensus

# Link in all the consensus files with samplename instead of lenghty bam samplename
for f in Projects/*/*.consensus.fasta
do
    # Base name of file
    bn=$(basename $f);
    # Link is relative from vcf_consensus folder to each consensus file
    ln -s ../${f} vcf_consensus/${bn%.bam.consensus.fasta}.fasta;
done
