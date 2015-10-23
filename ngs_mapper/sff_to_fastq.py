from ngs_mapper.reads import sffs_to_fastq
import os, sys
from functools import partial


def run(dir):
    sff_names = filter(lambda x: x.endswith('sff'), os.listdir(dir))
    sff_paths = map(partial(os.path.join, dir), sff_names)
    newext = lambda s: '.'.join(s.split('.')[:-1] + ['fastq'])
    outnames = map(newext, sff_paths)
    def wrapped_conv(a, b):
        print "Converting %s to fastq into %s" % (a, b)
        n = 0
        try:
            n=sffs_to_fastq([a], b, trim=True)
        except AssertionError as e: pass
        print "%s reads converted" % n
        return n
    return sum(map(wrapped_conv, sff_paths, outnames))

def main():
    dir = sys.argv[1]
    run(dir)

if __name__ == '__main__': main()



