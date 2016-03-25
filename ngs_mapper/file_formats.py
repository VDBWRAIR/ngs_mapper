from ngs_mapper.reads import sffs_to_fastq
import os, sys
from functools import partial
from glob import glob
import gzip
from Bio import SeqIO

{
    'sff' : 'sff',
    'ab1' : 'abi',  # SeqIO.convert(x, 'abi', y, 'fastq')
}

drop_ext = lambda s: '.'.join(s.split('.')[:-1])
swap_ext = lambda ext: lambda s:  drop_ext(s) + '.' + ext
find_ext = lambda ext: lambda dir: glob("%s/*%s" % (dir, ext))

def convert_sff(dir):
    sff_paths = find_ext('sff')(dir)
    outnames = map(swap_ext('fastq'), sff_paths)
    def wrapped_conv(a, b):
        print "Converting %s to fastq into %s" % (a, b)
        n = 0
        try:
            n=sffs_to_fastq([a], b, trim=True)
        except AssertionError as e:
            pass
        print "%s reads converted" % n
        return n
    return sum(map(wrapped_conv, sff_paths, outnames))

def convert_ab1(dir):
    for abi in find_ext('ab1')(dir):
        SeqIO.convert(abi, 'abi', swap_ext('fastq')(abi), 'fastq')

def convert_gzips(dir):
    for gz in find_ext('gz')(dir):
        dest = drop_ext(gz)
        with gzip.open( gz, 'rb' ) as input:
            with open(dest, 'w') as output:
                output.write(input.read())

def convert_formats(dir):
    convert_gzips(dir)
    convert_ab1(dir)
    convert_sff(dir)

def get_dir_arg():
    if os.path.isdir(sys.argv[1]):
        return sys.argv[1]
    else:
        raise ValueError("Path %s is not a directory" % sys.argv[1])

def main_convert_formats():
    convert_formats(get_dir_arg())


def main_sff_convert():
        convert_sff(get_dir_arg())

#        sff_names = filter(lambda x: x.endswith('sff'), os.listdir(dir))
#        sff_paths = map(partial(os.path.join, dir), sff_names)
#        run(sff_paths)





