include gmsl 
#TODO: add logging for bwa & !trimmomatic 
#TODO: fix basename calls?
#TODO: ngs_filter should really accept a list of files or individual files not a whole directory
# handle bwa keep-temp 
# this doesn't put stuff in tmpdir, This creates stuff in the directory it runs in ... 
#make will try to build any rule that's required even if it's LHS is empty

TRIMJAR := $(VIRTUAL_ENV)/lib/Trimmo*/*.jar
PAIR := *_R[12]_*.fastq
SFFS := $(wildcard *.sff)
FASTQS := $(filter-out unpaired.fastq, $(wildcard *.fastq) )
PAIRED :=  $(filter $(wildcard $(PAIR)), $(FASTQS)) 
FWD := $(wildcard filtered/*_R1_*.fastq)
REV := $(wildcard filtered/*_R2_*.fastq)

UNPAIRED :=  $(sort $(filter-out unpaired.fastq, $(filter-out $(wildcard $(PAIR)), $(FASTQS)) ))
#sort removes duplicates
F_UNPAIRED := $(sort $(addprefix filtered/filtered., $(UNPAIRED) $(SFFS:.sff=.fastq)))
F_PAIRED := $(sort $(addprefix filtered/filtered., $(PAIRED) ))
UNPTRIM := trimmed_reads/unpaired.trimmed.fastq
REF := {{ args.prefix }}.fasta 
#REF_INDEX := $(REF:=.amb) $(REF:=.ann) $(REF:=.bwt) $(REF:=.pac) $(REF:=.sa)
BAM := {{args.prefix}}.bam
BAI := $(BAM:.bam=.bai)

.PHONY:  all plots 
#TODO: make sure backup runs first
# backup always gets re-run, because F.fastq etc. go into the directory
# and they are alsways new
all:  plots flagstats.txt $(BAM).consensus.fasta

plots: {{args.prefix}}.reads.png  $(BAM).qualdepth.png $(BAM).qualdepth.json 
# this keeps getting re-run? 
# don't list directory as dependency
# prereqs might be newer
# unpaired.fastq keeps geting created should be in filtered 
#backup/backed: $(FASTQS) $(SFFS)
#	mkdir $(@D) && cp $^ $(@D) && touch $@

# Will it work if there is somename.sff + somename.fastq?

#
$(SFFS:.sff=.fastq): $(SFFS)
	mkdir backup && cp ./*.{fastq,sff} backup
	sff_to_fastq $^
	convert_read_quals $@ $(UNPAIRED)  #filter all unpaired
#track ngs_filter stats file

# http://www.cmcrossroads.com/article/trouble-hidden-targets
# http://stackoverflow.com/questions/7081284/gnu-make-multiple-targets-in-a-single-rule
# fixed using wildcard rule and specifying paired/unpaired within the action of the rule it depends on
filtered: $(PAIRED) $(UNPAIRED) $(SFFS:.sff=.fastq)
	@echo $(PAIRED) $(UNPAIRED)
	@ngs_filter . --outdir filtered {% if ngs_filter.dropNs %} --drop-ns  {% endif %} --platforms {{ ngs_filter.platforms|join(',') }} --index-min {{ ngs_filter.indexQualityMin }}

TMPFQ := filtered/up.tmp.fastq 
unpaired_trimmed: filtered
	mkdir -p trimmed_reads
	if [ "$(F_UNPAIRED)" ]; then cat $(F_UNPAIRED) > filtered/unpaired.fastq; java -jar $(TRIMJAR) SE -phred33 -trimlog trim.log filtered/unpaired.fastq $(TMPFQ) LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}};  fi;
	if [[ -s $(TMPFQ) ]]; then cp $(TMPFQ) $(UNPTRIM); fi;
	touch $@

trimargs_priv = $(1) $(1:=.unpaired) $(2) $(2:=.unpaired)  
trimouts =  $(call map, trimargs_priv, $(addprefix trimmed_reads/, $1 $2)) 
trim_paired = $(shell java -jar $(TRIMJAR) PE -trimlog trim.log $(1) $(2) $(call trimouts, $1, $2) LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}};) 

# $(wildcard) does not include arbitrary folders.
# '%' expansion is not lazy
# http://stackoverflow.com/questions/30303865/gnu-makes-prerequisite-expansion-is-not-lazy
trimmed_reads/pair: filtered
	mkdir -p trimmed_reads 
	@echo $(call pairmap, trim_paired, $(filter $(FWD), $(F_PAIRED)), $(filter $(REV), $(F_PAIRED))) && touch $@

$(BAI): $(BAM)
	samtools index $< 
	tagreads $< -CN {{ tagreads.CN }}

trimmed_reads/R.fastq:  trimmed_reads/pair
	cat trimmed_reads/*_R2_*.fastq > $@

trimmed_reads/F.fastq: trimmed_reads/pair
	cat trimmed_reads/*_R1_*.fastq > $@


#$(REF_INDEX): $(REF)
%.amb %.ann %.bwt %.pac %.sa: $(REF)
	bwa index $<

$(REF): {{ args.reference }}  
	@#if the reference is a directory, merge all fastas in it
	if [[ -d $< ]]; then cat $</* > $@; else cp $< $@; fi; 


PBAM := align/paired.bam
UNPBAM := align/unpaired.bam  #don't check ref_index anymore 
$(BAM): $(REF) trimmed_reads/F.fastq trimmed_reads/R.fastq unpaired_trimmed %.bwt  
	mkdir -p align
	bwa mem $(REF) $(word 2, $^) $(word 3, $^) > $(PBAM);
	if [[ -s $(UNPTRIM) ]];  then \
	   bwa mem $(REF) $(UNPTRIM) -t {{run_bwa_on_samplename.threads}} > $(UNPBAM);\
	   samtools merge - $(UNPBAM) $(PBAM) | samtools sort - $(basename $@);\
	else\
	   samtools sort $(PBAM) $(basename $@);\
	fi;

$(BAM).vcf: $(BAM) $(BAI) $(REF)
	base_caller $< $(word 3, $^) $@ -mith {{base_caller.minth}}

flagstats.txt: $(BAM) $(BAI)
	samtools flagstat $< > $@

$(BAM).qualdepth.png $(BAM).qualdepth.json: $(BAM) $(BAI) 
	graphsample $< -od $(@D) #directory of target 

{{args.prefix}}.reads.png:  unpaired_trimmed trimmed_reads/pair
	fqstats -o $@ trimmed_reads/*.fastq

$(BAM).consensus.fasta:   $(BAM).vcf
	vcf_consensus $< -i {{args.prefix}} -o $@ 
#
#TrimmomaticSE: Started with arguments: unpaired.fastq up.tmp.fastq LEADING:20 TRAILING:20 HEADCROP:0
#	Automatically using 8 threads
#Error: Unable to detect quality encoding
#	make: *** [unpaired_trimmed] Error 1
#make: Leaving directory `/home/AMED/michael.panciera/ngs_mapper-1.2/play/947'
#
#
# use directory instead of wildcards seems to work
# instead of trying to use wildcards could use github NGSXML strategy 
# of templating the rule for each individual file after ls'ing the directory.
# Converting 947__1__TI86__2012_08_06__Unk.sff to fastq into 947__1__TI86__2012_08_06__Unk.fastq
# /home/AMED/michael.panciera/.ngs/lib/python2.7/site-packages/Bio/SeqIO/SffIO.py:941: BiopythonParserWarning: Your SFF file is invalid, post index 4 byte null padding region contained data.
#   BiopythonParserWarning)
##   0 reads converted
#TrimmomaticPE: Started with arguments: -trimlog trim.log filtered/filtered.947_S32_L001_R1_001_2013_12_17.fastq filtered/filtered.947_S32_L001_R2_001_2013_12_17.fastq trimmed_reads/filtered/filtered.947_S32_L001_R1_001_2013_12_17.fastq trimmed_reads/filtered/filtered.947_S32_L001_R1_001_2013_12_17.fastq.unpaired trimmed_reads/filtered/filtered.947_S32_L001_R2_001_2013_12_17.fastq trimmed_reads/filtered/filtered.947_S32_L001_R2_001_2013_12_17.fastq.unpaired LEADING:20 TRAILING:20 HEADCROP:0
#TrimmomaticPE: Started with arguments: -trimlog trim.log 
#filtered/filtered.947_S32_L001_R1_001_2013_12_17.fastq 
#filtered/filtered.947_S32_L001_R2_001_2013_12_17.fastq 
#trimmed_reads/filtered/filtered.947_S32_L001_R1_001_2013_12_17.fastq 
#trimmed_reads/filtered/filtered.947_S32_L001_R1_001_2013_12_17.fastq.unpaired 
#trimmed_reads/filtered/filtered.947_S32_L001_R2_001_2013_12_17.fastq 
#trimmed_reads/filtered/filtered.947_S32_L001_R2_001_2013_12_17.fastq.unpaired 
#LEADING:20 TRAILING:20 HEADCROP:0
#Multiple cores found: Using 8 threads
#Quality encoding detected as phred33
#Exception in thread "main" java.io.FileNotFoundException: trimmed_reads/filtered/filtered.947_S32_L001_R1_001_2013_12_17.fastq (No such file or directory)
#        at java.io.FileOutputStream.open0(Native Method)
#        at java.io.FileOutputStream.open(FileOutputStream.java:270)
#        at java.io.FileOutputStream.<init>(FileOutputStream.java:213)
#        at java.io.FileOutputStream.<init>(FileOutputStream.java:162)
#        at org.usadellab.trimmomatic.fastq.FastqSerializer.open(FastqSerializer.java:28)
#        at org.usadellab.trimmomatic.TrimmomaticPE.process(TrimmomaticPE.java:275)
#        at org.usadellab.trimmomatic.TrimmomaticPE.run(TrimmomaticPE.java:498)
#        at org.usadellab.trimmomatic.Trimmomatic.main(Trimmomatic.java:35)
#
