cfg= --config ../config.yaml
trim_outdir = trimmed_reads
filtered_dir = fresh
samplename = 780
bamfile = $(addsuffix .bam, ${samplename})
vcf = $(addsuffix .vcf, ${bamfile})
consensus = $(addsuffix .consensus.fasta, ${bamfile})
flagstats = flagstats.txt
CN = ''
tdir = ${samplename}stats
head_crop = 0 
trim_qual = 20
bwalog = bwa.log
reference = 780/780.ref.fasta
minth = 0.8
${bamfile}: ${trim_outdir}
	run_bwa_on_samplename ${trim_outdir} ${reference} -o ${bamfile} ${cfg} > ${bwalog} 
tag_reads: ${bamfile}
	tagreads ${bamfile} -CN ${CN} ${cfg}

${vcf}:  tag_reads
	base_caller ${bamfile} ${reference} $@ -minth ${minth} ${cfg}

${flagstats}:
	samtools flagstat ${bamfile} > $@

${samplename}.bam.qualdepth.png ${samplename}.bam.qualdepth.json: ${vcf} ${flagstats}
	mkdir ${tdir}
	graphsample ${bamfile} -od ${tdir} 

$(basename ${bamfile}).reads.png: $(trim_outdir)/*.fastq
	echo trimming $<
	fqstats -o $@.reads.png $<

${samplename}.bam.consensus.fasta:  ${vcf} 
	vcf_consensus ${vcf} -i ${samplename} -o ${consensus} 
${trim_outdir} ${trim_outdir}/%.fastq: 
	trim_reads ${filtered_dir} -q ${trim_qual} -o ${trim_outdir} --head-crop ${head_crop} ${cfg} 

all: ${samplename}.bam.consensus.fasta info
	echo "donesies"

info: ${flagstats}  $(basename ${bamfile}).reads.png  ${samplename}.bam.qualdepth.png  ${samplename}.bam.qualdepth.json
