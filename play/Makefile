
PAIR := *_R[12]_*.fastq
SFFS = $(wildcard *.sff)
FASTQS := $(wildcard *.fastq)
UNPAIRED_ := $(filter-out $(wildcard $(PAIR)), $(FASTQS)) #this will get screwed up after sff conversion . . . 
PAIRED := $(filter $(wildcard $(PAIR)), $(FASTQS))
# this conj's the lists together
UNPAIRED = $(UNPAIRED_) $(SFFS:.sff=.fastq)

.PHONY: convert_sanger backup
all: backup trim 




all: $(PAIR)
	echo $+

unpaired.fastq: unpaired/*.fastq | con
	cat $+ > $@

unpaired/$(UNPAIRED): $(UNPAIRED)
	mkdir -p unpaired
	cp $+ unpaired 

	sffs_to_fastq unpaired/$+
	conv_read unpaired/$+

#trim_unpaired: trimmed/unpaired.fastq
#	if [ -s $< ]; then rm $<; fi;

convert_sff: unpaired/$(SFFS)
	if [ $+ ]; then sffs_to_fastq $+ ; fi; #current directory

trimmed/unpaired.fastq: unpaired.fastq
	mkdir -p trimmed
	java -jar trimmomatic.jar SE $< $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log


trimmed/$(PAIR): $(PAIR)
	mkdir -p trimmed
	java -jar trimmomatic.jar PE $+ $@ LEADING:{{trim_reads.q}} TRAILING:{{trim_reads.q}} HEADCROP:{{trim_reads.headcrop}} --trim_log out.log

backup: original/*
	mkdir -p original
	cp *.{fastq,sff} original


