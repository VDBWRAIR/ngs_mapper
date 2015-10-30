from Bio import SeqIO
import sys, os
from itertools import imap
import tempfile, shutil


def is_not_sanger(filepath):
    '''
    Inspect the top read in the file and see if the quality encoding
    max > 40
    '''
    #filepath = sys.argv[1]
    reads = SeqIO.parse( filepath, 'fastq' )
    r1 = next( reads )._per_letter_annotations['phred_quality']
    return not (max(r1) > 62)

def conv_read(_r):
    SANGER_OFFSET = 31
    conv = lambda x: ( x - SANGER_OFFSET )
    quals = _r._per_letter_annotations['phred_quality']
    fixed = map(conv, quals)
    _r._per_letter_annotations['phred_quality'] = fixed
    return _r

def convert_file(fn):
    reads = list(SeqIO.parse(fn, 'fastq'))
    os.remove(fn)
    with open(fn, 'w') as out:
        for rec in reads:
            out.write(rec.format("fastq-sanger"))

def convert_read_quals():
    files = sys.argv[1:]
    not_sangers = filter(is_not_sanger, files)
    for s in not_sangers:
        print "converting %s to sanger format." % s
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
