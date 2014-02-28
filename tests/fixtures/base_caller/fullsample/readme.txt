fullsample.sam is a subsampling from the 1090-01 sample so we have a bigger data set to test more thouroughly with

It was generated as follows:
samtools view -bhs 1.5 1090-01.bam > fullsample.bam

I think that would reduce the sample size by 50%. I then removed all of the unmapped reads in the sam file to futher reduce it's size.

I did some small tests with git and sam,bam and gzipped sam files and it
seemed that git was better at handling bam files and sam files so I decided to
do bam instead of sam for this instance even though in the long run if there
are a lot of updates to the bam file it may baloon the git repo size
