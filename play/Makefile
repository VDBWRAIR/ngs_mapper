#TODO: add logging for bwa & !trimmomatic 
# include tagreads
# handle bwa keep-temp 
# this doesn't put stuff in tmpdir
REF = {{ prefix }}.fasta
#aliasing a rule will cause stuff to be re-run even if it's .PHONY;the alternative is to store it in a variable!
REF_INDEX := $(REF).amb $(REF).ann $(REF).bwt $(REF).pac $(REF).sa 
PAIR := *_R[12]_*.fastq
SFFS = $(wildcard *.sff)
FASTQS := $(wildcard *.fastq)
UNPAIRED_ := $(filter-out unpaired.fastq, $(filter-out $(wildcard $(PAIR)), $(FASTQS))) #this will get screwed up after sff conversion . . . 
PAIRED := $(filter $(wildcard $(PAIR)), $(FASTQS)) 
BAM = {{ args.prefix }}.bam
BAI = $(BAM:.bam=.bai)

.PHONY:  align_unpaired
all: backup 

# This works even if there are no sff files (it doesn't get triggered?)
# Will it work if there is somename.sff + somename.fastq?
$(SFFS:.sff=.fastq): $(SFFS)
	@#mv $< $@
	sffs_to_fastq $^
	conv_read $(^:.sff=.fastq)

trimmed_reads/maybe.unpaired.fastq: $(UNPAIRED) $(SFFS:.sff=.fastq) 
	mkdir -p trimmed_reads
	if [ $^ ]; then  cat $^ > unpaired.fastq; rm $^;  java -jar trimmomatic.jar SE unpaired.fastq $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log; else touch $@; fi;

trimmed_reads/$(PAIR): $(PAIR)
	mkdir -p trimmed_reads
	java -jar trimmomatic.jar PE $^ $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log

backup: original/*
	mkdir -p original
	cp *.{fastq,sff} original 

#this shoudl work because align/paired.bam is listed first and should get run first
# $^ -> all prereqs
$(BAM): align/*.bam | align_unpaired
        samtools merge - $+ | samtools sort - $(basename $@)
	samtools sort $@ $(basename $@) 

$(BAI): $(BAM)
	samtools index $< 
	tagreads $< -CN {{ tagreads.CN }}

R.fastq:  trimmed_reads/*_R2_*.fastq
	cat $^ > $@

F.fastq: trimmed_reads/*_R1_*.fastq
	cat $^ > $@ 

$(REF_INDEX): $(REF)
	bwa index $<

$(REF): {{ args.reference }}  
	@#if the reference is a directory, merge all fastas in it
	if [[ -d $< ]]; then cat $</* > $@; else cp $< $@; fi;

#reference is first arg b/c third argument (reverse fastq) is optional
#This misses out on stuff like bwa_return_code which checks stdout/err to ensure that BWA ran correctly

bwa_mem = $(shell bwa mem $(2) $(3) $(1) -t {{run_bwa_on_samplename.threads}} > {{run_bwa_on_samplename.output}} ) 

$(BAM): F.fq P.fq trimmed_reads/maybe.unpairef.fastq $(REF) | unpaired.fastq
	# do all the merging and stuff here
	# don't let maybe.unpaired hang round in trimmed_read empty
align_unpaired: trimmed_reads/maybe.unpaired.fastq $(REF)     | indexed_ref
	if [ $< ]; then echo $(call bwa_mem, $(word 2, $^), $<); fi; 

align/paired.bam: $(REF) F.fastq R.fastq | indexed_ref
	echo $(call bwa_mem, $<, $(word 2, $^ ) $(word 3, $^))

$(BAM).vcf:  $(BAM) $(BAI) $(REF)
	base_caller $< $(word 3, $^) $@ -mith {base_caller.minth}

flagstats.txt: $(BAM) $(BAI)
	samtools flagstat $< > $@

$(BAM).qualdepth.png $(BAM).qualdepth.json: $(BAM) $(BAI) 
	graphsample $< -od $(@D) #directory of target

{{prefix}}.reads.png: trimmed_reads/*.fastq
	#ugh maybe.unpaired.fastq might be empty
	fqstats -o $@ $^

$(BAM).consensus.fasta:   $(BAM).vcf
	vcf_consensus $< -i {{args.prefix}} -o $@

