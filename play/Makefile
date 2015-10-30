include gsml

#TODO: add logging for bwa & !trimmomatic 

REF_INDEX := REF.amb REF.ann REF.bwt REF.pac REF.sa REF.amb 
PAIR := *_R[12]_*.fastq
SFFS = $(wildcard *.sff)
FASTQS := $(wildcard *.fastq)
UNPAIRED_ := $(filter-out $(wildcard $(PAIR)), $(FASTQS)) #this will get screwed up after sff conversion . . . 
PAIRED := $(filter $(wildcard $(PAIR)), $(FASTQS))
FORWARD := $(wildcard *_R1_*, $(PAIRED))
REVERSE := $(wildcard *_R2_*, $(PAIRED))
# watch out, gets applied to the list...

# this conj's the lists together
#UNPAIRED = $(filter-out unpaired.fastq, $(UNPAIRED_)) $(SFFS:.sff=.fastq)

#call is itself a funciton requires a comma after
# also just doing $(call func) w/in a rule will try to run the resul tof the call, which makes sense
UNPAIRED = $(filter-out unpaired.fastq, $(UNPAIRED_))
.PHONY:  align_unpaired
all: backup 

#make will try to build any rule that's required even if it's LHS is empty

# This works even if there are no sff files (it doesn't get triggered?)
# Will it work if there is somename.sff + somename.fastq?
$(SFFS:.sff=.fastq): $(SFFS)
	@#mv $< $@
	sffs_to_fastq $^
	conv_read $@

trimmed/maybe.unpaired.fastq: $(UNPAIRED) $(SFFS:.sff=.fastq) 
	mkdir -p trimmed_reads
	if [ $^ ]; then  cat $^ > unpaired.fastq; rm $^;  java -jar trimmomatic.jar SE unpaired.fastq $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log; else touch $@; fi;

trimargs_priv = $(1) $(1:=.unpaired) $(2) $(2:=.unpaired) 
trimouts =  $(call map, trimargs_priv, $(addprefix trimmed_reads/, $1 $2))
trim_paired = $(shell java -jar trimmomatic.jar PE $1 $2 $(call trimouts, $1, $2) LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log

trimmed_reads/$(PAIR): $(PAIR)
	mkdir -p trimmed_reads
	echo $(call pairmap, trim_paired, $^ ) 

backup: original/%
	mkdir -p original
	cp *.{fastq,sff} original 

#this shoudl work because align/paired.bam is listed first and should get run first
# $^ -> all prereqs
merged.bam: align/%.bam | align_unpaired
        samtools merge - $+ | samtools sort - $(basename $@)
	samtools sort $@ $(basename $@) 


# this works but just puts in the file name under RG name
#(.ngs)-bash-4.1$ samtools merge -rh rg.txt - s1.bam s.bam | samtools view -h - | head

merged.bai: merged.bam
	samtools index $< 

R.fastq:  trimmed/%_R2_%.fastq
	cat $^ > $@

F.fastq: trimmed/%_R1_%.fastq
	cat $^ > $@

#aliasing a rule will cause stuff to be re-run even if it's .PHONY; the alternative is to store it in a variable
$(REF_INDEX): REF.fasta
	bwa index $<

REF.fasta: {{ref}}
	if [[ -d $< ]]; then cat $</* > $@; else cp $< $@; fi;
#ordered-only will make Make ignore timestamp info . . . 
#this one is .PHONY
#reference is first arg b/c third argument (reverse fastq) is optional
bwa_mem = $(shell bwa mem $(2) $(3) $(1) -t {{bwa.threads}} > {{bwa.output}} )

# output file, keep_temp, threads 
# env = jinja2.Environment()
# env.from_string("-t {{bwa.threads}} -o").render( { "bwa" : { "threads" : 10 }})

align_unpaired:  REF.fasta trimmed/maybe.unpaired.fastq    | indexed_ref
	if [ $< ]; then bwa mem $< $(word $^ 2); fi;

#This misses out on stuff like bwa_return_code which ensure that BWA ran correctly

align/paired.bam: REF.fasta R.fastq P.fastq | indexed_ref
	$(call bwa_mem, $<, $(word $^ 2) $(word $^ 3))
# # can ngs_mapper take multilple sets of forward/reverse? yes
# do we really need to cat refs together? yes
