"""
Usage: make.py <readsdir> <reference> <prefix> [<target>]
"""
from docopt import docopt
from config import load_default_config
import sh
import sys
'''random utitlities'''
def flatten_dict(d, seperator='_'):
    def items():
        for key, value in d.items():
            try:
                for subkey, subvalue in flatten_dict(value, seperator=seperator).items():
                    if type(subvalue) is list: subvalue = ' '.join(map(str, subvalue))
                    yield key + seperator + subkey, subvalue
            except:
                yield key, value

    return dict(items())

def prepend(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

def strip_tag(s): return s[1:-1] if s.startswith('<') and s.endswith('>') else s
def stripargs(args_dict): return {strip_tag(k) : v for k, v in args_dict.items()}
def write(fname, content):
    with open(fname, 'w') as out:
        out.write(content)
    return fname
'''end random utilities'''

def run_make(makefile_str, target, filename='makefile.generated'):
    write(filename, makefile_str)
    log = open('pipeline.log', 'w')
    sh.make(target, f=filename,  _out=sys.stdout,  _err=log)
    log.close()


#MAKEFILE_TEMPLATE = 'Makefile.template'
make_template = '''
cfg= --config {CONFIG_PATH}
bwalog = bwa.log
all: consensus stats
consensus: %.consensus.fasta
alignment: %.bam
stats: {flagstats}  %.reads.png  %.qualdepth.png

%.bam: {trim_outdir}/%.fastq
        #bamfile = $(addsuffix .bam, {prefix})
	run_bwa_on_samplename {trim_outdir} {reference} -o $@ ${cfg} > ${bwalog}
	#TODO: tagreads should be a seperate step, even though it doesn't generate a new file
	tagreads $@ -CN {tagreads_CN_default} ${cfg}

%.vcf: %.bam
	base_caller $< {reference} $@ -minth {base_caller_minth_default} ${cfg}

{flagstats}:  %.bam
	samtools flagstat $< > $@

%.qualdepth.png %.qualdepth.json: %.bam %.vcf {flagstats}
	graphsample $< -od .

%.reads.png: $(trim_outdir)/*.fastq
	echo trimming $<
	fqstats -o $@  $<

%.consensus.fasta:  %.vcf
	vcf_consensus $< -i {prefix} -o $@

{trim_outdir}/%.fastq:  {readsdir}
	trim_reads $< -q {trim_reads_q_default} -o {trim_outdir} --head-crop {trim_reads_headcrop_default} ${cfg}
'''

def main():
    cfg_dict = load_default_config().yaml.yaml.yaml
    _vars = flatten_dict(cfg_dict, seperator='_')
    args = docopt(__doc__, version='0.1.0 alpha')
    cargs = stripargs(args)
    _vars.update(cargs)
    _vars['cfg'] = '{cfg}' #hacks
    _vars['bwalog'] = '{bwalog}' #more hacks
    _vars['CONFIG_PATH'] = '../config.yaml'
    _vars['trim_outdir'] = 'trimmed_dir'
    _vars['flagstats'] = 'flagstats.txt'
    makefile_str = make_template.format(**_vars)
    print makefile_str
    target = cargs['target'] or "all"
    run_make(makefile_str, target)

if __name__ == '__main__': main()

