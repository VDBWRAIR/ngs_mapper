from Bio import SeqIO
import sys, os
from itertools import imap 
import tempfile, shutil


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

def convert_file(fn):
    reads = list(SeqIO.parse(fn, 'fastq'))
    os.remove(fn)
    fixed = imap(conv_read, reads)
    SeqIO.write(fixed, fn, format='fastq')

def convert_sangers():
    files = sys.arv[1:]
    sangers = filter(is_sanger, files)
    for s in sangers:
        copy = create_temporary_copy(s)
        try:
            convert_file(s)
        except Exception, e:
            shutil.move(copy, s)
            sys.stderr.write( "Conversion of file %s failed with error %s" % (s, e)) 

def create_temporary_copy(path):
    '''http://stackoverflow.com/questions/6587516/how-to-concisely-create-a-temporary-file-that-is-a-copy-of-another-file-in-pytho'''
    temp_dir = tempfile.gettempdir()
    temp_path = os.path.join(temp_dir, 'temp_file_name')
    shutil.copy2(path, temp_path)
    return temp_path
