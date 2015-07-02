'''
Usage: ngs_filter <readdir> [--parallel]  [--drop-ns ] [--index-min=<index_min>] [--outdir <DIR>]

Options:
    --outdir=<DIR>,-o=<DIR>   outupt directory [Default: filtered]
    --parallel                Use python's multiprocessing to run on multiple cores

Help:
    --drop-ns                   Drop those reads which contain an N
    --index-min=<index_min>     Drop reads where the corresponding index is BELOW specified minimum.
'''
from functools import partial
import multiprocessing
from operator import methodcaller
from Bio import SeqIO
from docopt import docopt
from schema import Schema, Use, Optional
from itertools import ifilterfalse, izip
import re
import os
import sys
#TODO: should guarantee that there is an index file or allow there not to be one

''' given a fastq file and its index file (__R1__, __I1__), drop those reads that aren't good enough.'''
''' also, drop all reads with N value. '''
''' note: will have to deal with filename change. '''

get_index = partial(re.sub, r'_R([12])_', r'_I\1_')

def name_filtered(path, outdir):
    ''' rename with 'fitered.' prefix and inside the new path directory. '''
    rename = "filtered.{0}".format
    dirpath, filename = os.path.split(path)
    renamed = rename(filename)
    outd = outdir or '.'
    return os.path.join(outd, renamed)

def has_index(fn):
    ''' returns index path (ie _I1_ or _I2_)  or none.'''
    index = get_index(fn)
    return index if (os.path.exists(index) and index != fn) else None

def fqs_excluding_indices(readsdir):
    is_index = partial(re.search, r'_I[12]_')
    is_fastq = methodcaller('endswith', ('fq', 'fastq', 'sff'))
    def is_valid(fn):
        return is_fastq(fn) and not is_index(fn)
    all_files = os.listdir(readsdir)
    full_path = partial(os.path.join, readsdir)
    files = map(full_path, filter(is_valid, all_files))
    return files


def map_to_dir(readsdir, func, parallel=False):
    '''maps *func* to all fastq/sff files which are not indexes.
    fetch the fastqs and indexes of the directory and write the filtered results.'''
    files = fqs_excluding_indices(readsdir)
    if parallel:
        pool = multiprocessing.Pool()
        outpaths = pool.map(func, files)
        pool.close()
        pool.join()
        return outpaths
    else:
        return map(func, files)

def idx_filter(read, idxread, thresh):
    ''' AT or ABOVE threshold.'''
    return min(idxread._per_letter_annotations['phred_quality']) >= thresh

def make_filtered(readpath, idxQualMin, dropNs):
    ''' given a fastq file with an index, will filter on low-quality index entries, and drop all reads with N.
    If file does not have an index, only drops Ns.
    Raises an AssertionError if the index-quality filter drops all reads. '''
    index = has_index(readpath)
    if idxQualMin and not index:
        sys.stderr.write("Specified Non-null index quality minimum, but index for file {0} does not exist.\n".format(readpath))
    format = 'sff' if readpath.endswith('sff') else 'fastq'
    fq_open = partial(SeqIO.parse, format=format)
    if index and idxQualMin: 
         reads, idxreads = fq_open(readpath), fq_open(index)
         intermediate = (r for r, idx in izip(reads, idxreads) if idx_filter(r, idx, idxQualMin) )
    else:
        intermediate = fq_open(readpath)
    if dropNs:
        hasNs = lambda rec: 'N' in str(rec.seq).upper()
        return ifilterfalse(hasNs, intermediate)
    else:
        return intermediate

def write_filtered(readpath, idxQualMin, dropNs, outdir=None):
    '''write the results to the new directory'''
    results = make_filtered(readpath, idxQualMin, dropNs)
    outpath = name_filtered(readpath, outdir)
    results = list(results)
    sys.stderr.write('writing %s to %s' % (readpath, outpath)+'\n')
    sys.stderr.write(str(locals())+'\n')
    sys.stderr.write(open(readpath).read()[-10:] +'\n')


    num_written = SeqIO.write(results, outpath, 'fastq')
    assert num_written > 0, "Failed! Quality controls eliminated all reads. Drop-Ns was set to %s; \
        try again with lower quality min than %s. " %(dropNs, idxQualMin)
    return outpath

def write_post_filter(readsdir, idxQualMin, dropNs, outdir=None, parallel=False):
    '''execute write_filtered on the whole directory'''
    write_filters = partial(write_filtered, idxQualMin=idxQualMin, dropNs=dropNs, outdir=outdir)#, parallel=parallel)
    return map_to_dir(readsdir, write_filters, parallel)


def from_config(readsdir, config):
    ''' execute using the config dict derived from config.yaml '''
    print config
    qualmin = config['ngs_filter']['indexQualityMin']['default']
    dropNs = config['ngs_filter']['dropNs']['default']
    return write_post_filter(readsdir, qualmin, dropNs)

def mkdir_p(dir):
    ''' emulate bash command  $ mkdir -p '''
    if not os.path.exists(dir):
        os.mkdir(dir)

def main():
    scheme = Schema(
        { '<readdir>' : os.path.isdir,
         Optional('--drop-ns') : bool,
         Optional('--parallel') : bool,
         Optional('--index-min') : Use(int, error="--index-min takes integer"),
         '--outdir' : str
         })

    raw_args = docopt(__doc__, version='Version 0')
    args = scheme.validate(raw_args)
    mkdir_p(args['--outdir'])
    dropNs, idxMin = args['--drop-ns'], args['--index-min']
    if not (dropNs or idxMin > -1):
        raise ValueError("No filter specified, drop Ns:%s, Index Quality Min:%s" % (dropNs, idxMin))
    minmin, minmax = -1, 50
    if not ( minmin <= idxMin <= minmax):
        raise ValueError("Invalid Index Quality Minimum specified: %s  is not a valid value between %s and %s" (idxMin, minmin, minmax))
    sys.stderr.write( "\nfiltering with specifications, drop Ns:%s, Index Quality Min:%s\nfrom folder %s to folder %s\n" % (dropNs, idxMin, args['<readdir>'], args['--outdir']))
    sys.stderr.write(os.getcwd())
    sys.stderr.write('\nparallel: %s\n' % args['--parallel'])
    outpaths = write_post_filter(args['<readdir>'], idxMin, dropNs, args['--outdir'], args['--parallel'])
    return 0


