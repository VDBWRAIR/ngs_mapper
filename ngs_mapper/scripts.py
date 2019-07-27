from Bio import SeqIO
import sys, os
from itertools import imap
def is_sanger():
    '''
    Inspect the top read in the file and see if the quality encoding
    max > 40
    '''
    filepath = sys.argv[1]
    reads = SeqIO.parse( filepath, 'fastq' )
    r1 = next( reads )._per_letter_annotations['phred_quality']
    return max(r1) > 40

def conv_read(r):
    SANGER_OFFSET = 33
    conv = lambda x: chr( int(x) - SANGER_OFFSET )
    quals = r._per_letter_annotations['phred_quality']
    fixed = map(conv, quals)
    r._per_letter_annotations['phred_quality'] = fixed

def convert_sangers():
    files = sys.arv[1:]
    sangers = filter(is_sanger, files)
    def convert_file(fn):
        reads = list(SeqIO.parse(fn, 'fastq'))
        os.remove(fn)
        fixed = imap(conv_read, reads)
        SeqIO.write(fixed, fn, format='fastq')
    for s in sangers:
        convert_file(s)
