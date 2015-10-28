#TODO: add logging for bwa & !trimmomatic 
# include tagreads
# handle bwa keep-temp

REF_INDEX := REF.amb REF.ann REF.bwt REF.pac REF.sa REF.amb 
PAIR := *_R[12]_*.fastq
SFFS = $(wildcard *.sff)
FASTQS := $(wildcard *.fastq)
UNPAIRED_ := $(filter-out unpaired.fastq, $(filter-out $(wildcard $(PAIR)), $(FASTQS))) #this will get screwed up after sff conversion . . . 
PAIRED := $(filter $(wildcard $(PAIR)), $(FASTQS)) 

.PHONY:  align_unpaired
all: backup 

# This works even if there are no sff files (it doesn't get triggered?)
# Will it work if there is somename.sff + somename.fastq?
$(SFFS:.sff=.fastq): $(SFFS)
	@#mv $< $@
	sffs_to_fastq $^
	conv_read $(^:.sff=.fastq)

trimmed/maybe.unpaired.fastq: $(UNPAIRED) $(SFFS:.sff=.fastq) 
	mkdir -p trimmed
	if [ $^ ]; then  cat $^ > unpaired.fastq; rm $^;  java -jar trimmomatic.jar SE unpaired.fastq $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log; else touch $@; fi;

trimmed/$(PAIR): $(PAIR)
	mkdir -p trimmed
	java -jar trimmomatic.jar PE $^ $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log

backup: original/%
	mkdir -p original
	cp *.{fastq,sff} original 

#this shoudl work because align/paired.bam is listed first and should get run first
# $^ -> all prereqs
merged.bam merged.bai: align/%.bam | align_unpaired
        samtools merge - $+ | samtools sort - $(basename $@)
	samtools sort $@ $(basename $@) 

merged.bai: merged.bam
	samtools index $< 
	tagreads $< -CN {{ tagreads.CN }}

R.fastq:  trimmed/%_R2_%.fastq
	cat $^ > $@

F.fastq: trimmed/%_R1_%.fastq
	cat $^ > $@ 

#aliasing a rule will cause stuff to be re-run even if it's .PHONY;the alternative is to store it in a variable!
$(REF_INDEX): REF.fasta
	bwa index $<

REF.fasta: {{ref}}  
	@#if the reference is a directory, merge all fastas in it
	if [[ -d $< ]]; then cat $</* > $@; else cp $< $@; fi;

#reference is first arg b/c third argument (reverse fastq) is optional
#This misses out on stuff like bwa_return_code which ensure that BWA ran correctly
bwa_mem = $(shell bwa mem $(2) $(3) $(1) -t {{bwa.threads}} > {{bwa.output}} ) 

align_unpaired:  REF.fasta trimmed/maybe.unpaired.fastq    | indexed_ref
	if [ $< ]; then echo $(call bwa_mem, $<, $(word 2, $^)); fi; 

align/paired.bam: REF.fasta R.fastq P.fastq | indexed_ref
	echo $(call bwa_mem, $<, $(word 2, $^ ) $(word 3, $^))
