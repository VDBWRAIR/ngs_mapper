import sys
import os
from os.path import *
import subprocess
import argparse

import bqd, graph_qualdepth as qd
from bam_to_qualdepth import set_unmapped_mapped_reads
import json
import log

logc = log.get_config( 'graphsample.log' )
logger = log.setup_logger( 'graphsample', logc )

def main():
    args = parse_args()
    args = handle_args( args )
    if not args.qualdepth:
        jfile = make_json( args.bamfile, args.outpath )
    else:
        jfile = args.qualdepth
    pngfile = make_image( jfile, args.outpath )

def make_json( bamfile, outpathprefix ):
    pileup = bqd.mpileup( bamfile )
    stats = bqd.parse_pileup( pileup )
    set_unmapped_mapped_reads( bamfile, stats )
    outfile = outpathprefix + '.qualdepth.json'
    with open( outfile, 'w' ) as fh:
        json.dump( stats, fh )

    return outfile

def make_image( jfile, outpathprefix ):
    prefix = basename( outpathprefix )
    imgdir = join( dirname(outpathprefix), 'qualdepth' )
    if not exists( imgdir ):
        os.mkdir( imgdir )
    outfile = join( imgdir, basename(outpathprefix) + '.qualdepth.' )
    j = json.load( open(jfile) )
    imagelist = []
    for ref in [r for r in j if r != 'unmapped_reads']:
        refname=normalize_ref(ref)
        title = prefix + ' ' + refname
        of = outfile + refname + '.png'
        qd.make_graphic( jfile, of, ref=ref, titleprefix=title )
        imagelist.append( of )
    imagelist.append( outpathprefix + '.qualdepth.png' )
    run_montage( *imagelist, compress='JPEG', quality=25, geometry='+1+1' )
    return outfile

def normalize_ref( refname ):
    '''
        Replace punctuation with _
    '''
    import string
    name = ''
    for c in refname:
        if c in string.punctuation + string.whitespace:
            name += '_'
        else:
            name += c
    return name

def run_montage( *args, **kwargs ):
    '''
        Runs montage on all the images to create the
        final qualdepth.png file
        
        @param outprefix - Prefix that is common to all files to montage together

        @param args - Files to montage with the last item in the list being the output file
        @param kwargs - options to pass to montage(compress='JPEG',quality=25)

        @returns the output file path which should be outprefix + '.png'
    '''
    cmd = ['montage']
    for k,v in kwargs.items():
        cmd += ['-{}'.format(k),str(v)]
    cmd += args
    logger.debug( 'Running {}'.format( ' '.join(cmd) ) )
    try:
        subprocess.check_call( cmd )
    except subprocess.CalledProcessError as e:
        logger.critical('Could not build montage image because {0}'.format(e))
    except OSError as e:
        logger.critical('Montage command is missing')
    return args[-1]
    
def handle_args( args ):
    if args.outprefix is not None:
        outprefix = args.outprefix
    else:
        bamname = basename( args.bamfile )
        outprefix = bamname

    args.outpath = join( args.outdir, outprefix )

    return args

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Runs bam_to_qualdepth as well as graph_qualdepth from BamCoverage'
    )

    parser.add_argument(
        'bamfile',
        help='Path to bamfile'
    )

    parser.add_argument(
        '-od',
        default=os.getcwd(),
        dest='outdir',
        help='Where to place the output files(outdir/outprefix.{png,json})[Default: Current directory]'
    )

    parser.add_argument(
        '-op',
        '--out-prefix',
        dest='outprefix',
        default=None,
        help='How to name the output files. Default is to prefix with the bamfile name'
    )

    parser.add_argument(
        '-qualdepth',
        dest='qualdepth',
        default=None,
        help='Specify an already existing qualdepth.json file so it doesn\'t have to be recreated'
    )

    return parser.parse_args( args )
