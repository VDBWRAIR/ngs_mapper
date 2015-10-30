

concat lists by simply placing them together
LISTS = $(L1) $(L2)

note that the commas are weird (but consistent) with function calls:
gsml includes map, pairedmap (zip-with) etc.

$(addprefix <nocomma> arg1 <commas...> arg2)
$(call func <nocomma> arg1, arg2)
$(shell <no-commas-ever)
#call is itself a funciton requires a comma after
$(var:=.newsuff)  -- add a suffix
$(var:.old=.new) -- replace a suffix
=======

#aliasing a rule will cause stuff to be re-run even if it's .PHONY;the alternative is to store it in a variable!
reference is first arg b/c third argument (reverse fastq) is optional
 samtools sort creates a new file doesn't touch existing file
This misses out on stuff like bwa_return_code which checks stdout/err to ensure that BWA ran correctly
first arg is output lol
 
 In a pattern rule that has multiple targets (see Introduction to Pattern Rules), ‘$@’ is the name of whichever target caused the rule’s recipe to be run. 
ordered-only will make Make ignore timestamp info . . . 
 this conj's the lists together
UNPAIRED = $(filter-out unpaired.fastq, $(UNPAIRED_)) $(SFFS:.sff=.fastq)
aliasing a rule will cause stuff to be re-run even if it's .PHONY;the alternative is to store it in a variable!
`jinja2.Environment().from_string(s).render(dict)`
make will try to build any rule that's required even if it's LHS is empty
call is itself a funciton requires a comma after
also just doing $(call func) w/in a rule will try to run the resul tof the call, which makes sense

this below works but just puts in the file name under RG name
1$ samtools merge -rh rg.txt - s1.bam s.bam | samtools view -h - | head
= is lazy, := is strict
 when to remove empty unpaired file?
problems: inconvenient to have a dependency which is not used in the rule because of how $+ works
#negative regex for unpaired
#variables (required for negative regex) are instantiated at startup?
#
 | -> order-only pre-requisite, useful for converting the unpaired file

note: if a pre-req matches multiple files, you still need to use $+ to grab all of it!
sff->fastq
run ngs_filter 

lots of file name-fu
%.bam: unpaired.fastq paried.fastq
     if [-s unpaired.fastq] run_bwa   $<
     run_bwa $>

merged.bam: *.bam
    samtools merge -r -h rg.txt $+ > $@

 $* automatic variable, which is defined as the part of the filename that matched the % character in a pattern rule.

$*_R1_.fastq
http://stackoverflow.com/questions/5561123/how-to-get-a-particular-dependency-from-the-dependency-list-for-a-target?lq=1

make -C subdir/Makefile will cd to subdir first then run make
---trim reads
-capture trim stats
-SE/PE option
-detect phred33 
merge non-empty unpaired files

phred33 can be converted ez
convert any phred33 files
merge unpaired
think of little tools like little bash tools 

trim_reads
================

don't trim index reads, only trim fastq/sff

(used to filter platforms)
trim the reads

merge non-empty unpaired files into one file

create trim and stats dirs

handle the SE/PE option in trimmomatic
run trimmomatic
write stats file

detect phred33 if data.is_sanger_readfile(args[1]): and add -phred33 if it's there

PE is paired-end, SE is single

create stepname:value string, get jar path and run
apparently don't run cutadapt

opts:
headcrop,  -q (qual threshold)


run_bwa_on_sample_name
===========================
{:-t threads :reference :reads :platforms?

compile_reads
Compiles all read files inside of readfilelist into respective files.
    Creates F.fq, R.fq and/or NP.fq depending on the reads found in readfilelist
    Only compiles fastq files. If others are given an exception will be raised

 compiles_refs -> just cats refs together if "ref is  a dir"

 run bwa: index the reference, then run bwa

run bwa with paired reads if exists, unpaired read if exists (done separately!)  -> merge bams

index the bam file


tagreads
========
[:bamfiles :SM :CN]
 SM is 
can be replaced by:

samtools merge -r -h? 
perl -e 'print "@RG\tID:ga\tSM:hs\tLB:ga\tPL:Illumina\n@RG\tID:454\tSM:hs\tLB:454\tPL:454\n"' > rg.txt
samtools merge -rh rg.txt - ga.bam 454.bam | samtools rmdup - - | samtools rmdup -s - aln.bam
seems to make a new bam, sort and replace it


base_caller
===========
leave alone



 lots of file name-fu
can ngs_mapper take multilple sets of forward/reverse? yes
 do we need to cat refs together? yes
