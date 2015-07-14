'''
Usage: ngs_filter <readdir> [--parallel]  [--drop-ns ] [--index-min=<index_min>] [--outdir <DIR>]

Options:
    --outdir=<DIR>,-o=<DIR>   outupt directory [Default: filtered]
    --parallel                Use python's multiprocessing to run on multiple cores

Help:
    If an argument is not given for a parameter, that filter is not applied. If no filter parameters are provided, an error is raised.
    --drop-ns                   Drop those reads which contain an N
    --index-min=<index_min>     Drop reads where the corresponding index is BELOW specified minimum. Must be between 1 and 50.
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
import warnings
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
    renamed = renamed[:-3] + 'fastq' if renamed.endswith('sff') else renamed
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
    if not files:
        raise ValueError("No fastq or sff files found in directory %s" % readsdir)
    if parallel:
        pool = multiprocessing.Pool()
        outpaths = pool.map(func, files)
        pool.close()
        pool.join()
        return outpaths
    else:
        print "mapping filters over read files %s in directory %s" % (files, readsdir)
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
    try:
        num_written = SeqIO.write(results, outpath, 'fastq')
        print "filtered reads from %s will be written to %s" % (readpath, outpath)
        print "%s reads left after filtering." % num_written
        if  num_written <= 0:
            warnings.warn("No reads left after filtering! Quality controls eliminated all reads. Drop-Ns was set to %s; maybe try again with lower quality min than %s. " %(dropNs, idxQualMin))
    except AssertionError, E:
        sys.stderr.write(str(E))
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
         Optional('--index-min') : Use(lambda x: int(x) if x else x, error="--index-min expects an integer"),
         '--outdir' : str
         })

    raw_args = docopt(__doc__, version='Version 0')
    args = scheme.validate(raw_args)
    mkdir_p(args['--outdir'])
    dropNs, idxMin = args['--drop-ns'], args['--index-min']
    minmin, minmax = 1, 50
    if not (dropNs or (minmin <= idxMin <=minmax)):
        raise ValueError("No filter specified, drop Ns:%s, Index Quality Min:%s" % (dropNs, idxMin))
    status = "\nfiltering with specifications, drop Ns:%s, Index Quality Min:%s\nfrom folder %s to folder %s" % (dropNs, idxMin, args['<readdir>'], args['--outdir'])
    print status
    outpaths = write_post_filter(args['<readdir>'], idxMin, dropNs, args['--outdir'], args['--parallel'])
    return 0


