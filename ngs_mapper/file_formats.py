from ngs_mapper.reads import sffs_to_fastq
import os, sys
from functools import partial
from glob import glob
import gzip
from Bio import SeqIO
from . import log

logger = log.setup_logger(__name__, log.get_config())

{
    'sff' : 'sff',
    'ab1' : 'abi',  # SeqIO.convert(x, 'abi', y, 'fastq')
}

drop_ext = lambda s: '.'.join(s.split('.')[:-1])
swap_ext = lambda ext: lambda s:  drop_ext(s) + '.' + ext
find_ext = lambda ext: lambda dir: glob("%s/*%s" % (dir, ext))
swap_dir = lambda dir: lambda p: os.path.join(dir, basename(p))

def convert_sff(dir, outdir):
    sff_paths = find_ext('sff')(dir)
    outnames = map(swap_ext('fastq'), sff_paths)
    outnames = map(swap_dir(outdir), outnames)
    def wrapped_conv(a, b):
        logger.info('Converting {0} to {1}'.format(a, b))
        n = 0
        try:
            n=sffs_to_fastq([a], b, trim=True)
        except AssertionError as e:
            pass
        logger.info("{0} reads converted".format(n))
        return n
    return sum(map(wrapped_conv, sff_paths, outnames))

def convert_ab1(dir, outdir):
    for abi in find_ext('ab1')(dir):
        dest = swap_ext('fastq')(abi)
        dest = swap_dir('outdir')(dest)
        logger.info('Converting {0} to {1}'.format(abi, dest))
        SeqIO.convert(abi, 'abi', dest, 'fastq')

def convert_gzips(dir, outdir):
    for gz in find_ext('gz')(dir):
        dest = swap_dir(outdir)(drop_ext(gz))
        with gzip.open( gz, 'rb' ) as input:
            with open(dest, 'w') as output:
                logger.info('Unpacking {0} to {1}'.format(gz, dest))
                output.write(input.read())

def convert_formats(dir, outdir):
    convert_gzips(dir, outdir)
    convert_ab1(dir, outdir)
    convert_sff(dir, outdir)

def get_dir_args():
    if os.path.isdir(sys.argv[1]) and  os.path.isdir(sys.argv[2]):
        return sys.argv[1], sys.argv[2]
    else:
        raise ValueError("Path %s or %s is not a directory" % (sys.argv[1], sys.argv[2]))

def main_convert_formats():
    convert_formats(*get_dir_args())


def main_sff_convert():
        convert_sff(*get_dir_args())

#        sff_names = filter(lambda x: x.endswith('sff'), os.listdir(dir))
#        sff_paths = map(partial(os.path.join, dir), sff_names)
#        run(sff_paths)





