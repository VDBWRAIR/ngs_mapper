include gsml 
#TODO: add logging for bwa & !trimmomatic 
#TODO: fix basename calls?
# handle bwa keep-temp 
# this doesn't put stuff in tmpdir, This creates stuff in the directory it runs in ... 
#make will try to build any rule that's required even if it's LHS is empty

TRIMJAR := $(VIRTUAL_ENV)/lib/Trimmo*/*.jar
PAIR := *_R[12]_*.fastq
SFFS := $(wildcard *.sff)
FASTQS := $(wildcard *.fastq) 
PAIRED :=  $(filter $(wildcard $(PAIR)), $(FASTQS)) 
FORWARD := $(wildcard *_R1_*, $(PAIRED))
REVERSE := $(wildcard *_R2_*, $(PAIRED))
UNPAIRED :=  $(filter-out $(wildcard $(PAIR)), $(FASTQS)) 
#sort removes duplicates
F_UNPAIRED := $(sort $(addprefix filtered/filtered., $(UNPAIRED) $(SFFS:.sff=.fastq)))
F_PAIRED := $(sort $(addprefix filtered/filtered., $(PAIRED) $(SFFS:.sff=.fastq)))
UNPTRIM := trimmed_reads/unpaired.trimmed.fastq
REF := {{ args.prefix }}.fasta 
REF_INDEX := $(REF:=.amb) $(REF:=.ann) $(REF:=.bwt) $(REF:=.pac) $(REF:=.sa)
BAM := {{args.prefix}}.bam
BAI := $(BAM:.bam=.bai)

do: $(F_UNPAIRED) filtered/$(PAIR)
	echo $(F_UNPAIRED) filtered/$(PAIR)

.PHONY:  all plots 
#TODO: make sure backup runs first
all: backup plots flagstats.txt $(BAM).consensus.fasta 

plots: {{args.prefix}}.reads.png  $(BAM).qualdepth.png $(BAM).qualdepth.json 

backup: $(FASTQS) $(SFFS)
	mkdir $@ && cp $^ $@

# This works even if there are no sff files (it doesn't get triggered?)
# Will it work if there is somename.sff + somename.fastq?
$(SFFS:.sff=.fastq): $(SFFS)
	sffs_to_fastq $^
	conv_read $@ $(UNPAIRED) 

trimmed/maybe.unpaired.fastq: $(UNPAIRED) $(SFFS:.sff=.fastq) 
	mkdir -p trimmed_reads
	if [ $^ ]; then  cat $^ > unpaired.fastq; rm $^;  java -jar trimmomatic.jar SE unpaired.fastq $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log; else touch $@; fi;

$(F_UNPAIRED) $(F_PAIRED): $(PAIRED) $(UNPAIRED)
	ngs_filter . --outdir $(@D) {% if ngs_filter.dropNs %} --drop-ns  {% endif %} --platforms {{ ngs_filter.platforms|join(',') }} --index-min {{ ngs_filter.indexQualityMin }}

TMPFQ := up.tmp.fastq 
unpaired_trimmed: $(F_UNPAIRED) 
	mkdir -p trimmed_reads
	if [ "$^" ]; then cat $^ > unpaired.fastq; rm $^;  java -jar $(TRIMJAR) SE unpaired.fastq $(TMPFQ) LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log;  fi;
	if [[ -s $(TMPFQ) ]]; then cp $(TMPFQ) $(UNPTRIM); fi;
	touch $@

trimargs_priv = $(1) $(1:=.unpaired) $(2) $(2:=.unpaired) 
trimouts =  $(call map, trimargs_priv, $(addprefix trimmed_reads/, $1 $2))
trim_paired = $(shell java -jar $(TRIMJAR) PE $1 $2 $(call trimouts, $1, $2) LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log trim.log)
trimmed_reads/$(PAIR): $(F_PAIRED)
	mkdir -p trimmed_reads
	echo $(call pairmap, trim_paired, $^ ) 

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


PBAM := align/paired.bam
UNPBAM := align/unpaired.bam 
$(BAM): $(REF) F.fq P.fq unpaired_trimmed/ $(REF_INDEX)
	mkdir -p align
	bwa mem $(REF) $(word 2, $^) $(word 3 $^) > $(PBAM);
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

{{args.prefix}}.reads.png:  unpaired_trimmed trimmed_reads/$(PAIR) 
	fqstats -o $@ trimmed_reads/*.fastq

$(BAM).consensus.fasta:   $(BAM).vcf
	vcf_consensus $< -i {{args.prefix}} -o $@ 
