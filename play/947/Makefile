#TODO: add logging for bwa & !trimmomatic 
# handle bwa keep-temp 
# this doesn't put stuff in tmpdir
# This creates stuff in the directory it runs in ... 
TRIM_JAR := $(VIRTUAL_ENV)/lib/Trimmo*/*.jar
PAIR := *_R[12]_*.fastq
SFFS := $(wildcard *.sff)
FASTQS := $(wildcard *.fastq)
PAIRED :=  $(filter $(wildcard $(PAIR)), $(FASTQS)) 
#sort removes duplicates
UNPAIRED := $(sort $(filter-out $(wildcard $(PAIR)), $(FASTQS)) )
F_UNPAIRED := $(addprefix filtered/filtered., $(UNPAIRED) $(SFFS:.sff=.fastq))
F_PAIRED := $(sort $(addprefix filtered/filtered., $(PAIRED) $(SFFS:.sff=.fastq)))
UNPTRIM := trimmed_reads/unpaired.trimmed.fastq
REF := 947.fasta 
REF_INDEX := $(REF:=.amb) $(REF:=.ann) $(REF:=.bwt) $(REF:=.pac) $(REF:=.sa)
BAM := 947.bam
BAI := $(BAM:.bam=.bai)


.PHONY:  all plots 
#make sure backup runs first
all: backup plots flagstats.txt $(BAM).consensus.fasta

plots: 947.reads.png  $(BAM).qualdepth.png $(BAM).qualdepth.json 
#TODO: fix basename calls?
backup: $(FASTQS) $(SFFS)
	mkdir $@ && cp $^ $@

# This works even if there are no sff files (it doesn't get triggered?)
# Will it work if there is somename.sff + somename.fastq?
$(SFFS:.sff=.fastq): $(SFFS)
	@#mv $< $@
	sffs_to_fastq $^
	conv_read $@ $(UNPAIRED) 

$(F_UNPAIRED) $(F_PAIRED): $(PAIRED) $(UNPAIRED)
	ngs_filter . --outdir $(@D)  --platforms MiSeq,Sanger,Roche454,IonTorrent --index-min 1


TMPFQ := up.tmp.fastq 
unpaired_trimmed: $(F_UNPAIRED) 
	mkdir -p trimmed_reads
	if [ "$^" ]; then cat $^ > unpaired.fastq; rm $^;  java -jar $(TRIM_JAR) SE unpaired.fastq $(TMPFQ) LEADING:20 TRAILING:20 HEADCROP:0 --trim_log out.log;  fi;
	if [[ -s $(TMPFQ) ]]; then cp $(TMPFQ) $(UNPTRIM); fi;
	touch $@

trimmed_reads/$(PAIR): $(F_PAIRED)
	@echo "trimming paired reads"
	mkdir -p $(@D)
	java -jar $(TRIM_JAR) PE $^ $(addprefix trimmed_reads/, $^)  LEADING:20 TRAILING:20 HEADCROP:0 --trim_log out.log

$(BAI): $(BAM)
	samtools index $< 
	tagreads $< -CN None

R.fastq:  trimmed_reads/*_R2_*.fastq
	cat $^ > $@

F.fastq: trimmed_reads/*_R1_*.fastq
	cat $^ > $@ 

$(REF_INDEX): $(REF)
	bwa index $<

$(REF): 947.ref.fasta  
	@#if the reference is a directory, merge all fastas in it
	if [[ -d $< ]]; then cat $</* > $@; else cp $< $@; fi;


PBAM := align/paired.bam
UNPBAM := align/unpaired.bam 
$(BAM): $(REF) F.fq P.fq unpaired_trimmed/ $(REF_INDEX)
	mkdir -p align
	bwa mem $(REF) $(word 2, $^) $(word 3 $^) > $(PBAM);
	if [[ -s $(UNPTRIM) ]];  then \
	   bwa mem $(REF) $(UNPTRIM) -t 1 > $(UNPBAM);\
	   samtools merge - $(UNPBAM) $(PBAM) | samtools sort - $(basename $@);\
	else\
	   samtools sort $(PBAM) $(basename $@);\
	fi;

$(BAM).vcf: $(BAM) $(BAI) $(REF)
	base_caller $< $(word 3, $^) $@ -mith 0.8

flagstats.txt: $(BAM) $(BAI)
	samtools flagstat $< > $@

$(BAM).qualdepth.png $(BAM).qualdepth.json: $(BAM) $(BAI) 
	graphsample $< -od $(@D) #directory of target 

947.reads.png:  unpaired_trimmed trimmed_reads/$(PAIR) 
	fqstats -o $@ trimmed_reads/*.fastq

$(BAM).consensus.fasta:   $(BAM).vcf
	vcf_consensus $< -i 947 -o $@ 
