"""
Here be some documentation
"""

import subprocess
import os
import argparse
import sys
from os.path import basename, join, isdir, dirname, expandvars
from glob import glob
import tempfile
import reads
import shlex
import data
from ngs_mapper import compat

import log
lconfig = log.get_config()
logger = log.setup_logger( 'trim_reads', lconfig )

def main():
    args = parse_args()
    trim_reads_in_dir(
        args.readsdir,
        args.q,
        args.outputdir,
        head_crop=args.headcrop,
        platforms=args.platforms,
        primer_info=[args.primer_file, args.primer_seed, args.palindrom_clip, args.simple_clip]
    )

def trim_reads_in_dir( *args, **kwargs ):
    '''
        Trims all read files in a given directory and places the resulting files into out_path directory
        Combines all unpaired trimmed fastq files into a single file as a result of running Trimmomatic

        :param str readdir: Directory with read files in it(sff and fastq only)
        :param int qual_th: What to pass to cutadapt -q
        :param str out_path: Output directory path
        :param int head_crop: How many bases to crop off ends
        :param list platforms: List of platform's reads to use
    '''
    readdir = args[0]
    qual_th = args[1]
    out_path = args[2]
    headcrop = kwargs.get('head_crop', 0)
    platforms = kwargs.get('platforms', None)
    primer_info = kwargs.get('primer_info')
    logger.info(
        "Only accepting the following platform's read files: {0}".format(
            platforms
        )
    )

    # Only sff and fastq files
    #reads = [f for f in os.listdir(readdir) if f.endswith('sff') or f.endswith('fastq')]
    platreads = data.reads_by_plat( readdir )
    # Make out_path
    if not isdir( out_path ):
        os.mkdir( out_path )
    # Trim all the reads
    unpaired = []
    for plat, reads in platreads.iteritems():
        if platforms is None or plat in platforms:
            for r in reads:
                if '_I1_' in r or '_I2_' in r:
                    logger.debug("Skipped MiSeq index read {0}".format(r))
                    continue
                inreads = None
                if isinstance(r,str):
                    # Only accept .sff and .fastq
                    if r.endswith('.sff') or r.endswith('.fastq'):
                        inreads = r
                        outreads = join(out_path, basename(r).replace('.sff','.fastq'))
                else:
                    inreads = r
                    outreads = [join(out_path, basename(pr)) for pr in r]
                # Sometimes inreads not set because .ab1 files
                if inreads is None:
                    continue
                try:
                    r = trim_read( inreads, qual_th, outreads, head_crop=headcrop, primer_info=primer_info )
                    logger.debug("Output from trim_read {0}".format(r))
                    unpaired += r[1::2]
                    logger.debug("Added {0} to unpaired list".format(r[1::2]))
                except subprocess.CalledProcessError as e:
                    print e.output
                    raise e
        else:
            logger.info("{0} are being excluded as they are not in {1}".format(
                reads,platforms
            ))
    #!!
    # Combine all *.fastq.unpaired into one file for mapping as SE
    #!!
    # Filter out unpaired files that are empty
    notempty = filter( lambda f: os.stat(f).st_size > 0, unpaired )
    if notempty:
        logger.info("Combining all unpaired trimmed files into a single file")
        out_unpaired = join( out_path, 'unpaired_trimmed.fastq' )
        with open( out_unpaired, 'w' ) as fw:
            for up in notempty:
                with open(up) as fr:
                    fw.write( fr.read() )
    else:
        logger.debug("All unpaired trimmed files are empty")
    # Remove the read now as it is no longer needed
    for up in unpaired:
        os.unlink( up )

def trim_read( *args, **kwargs ):
    '''
        Trims the given readpath file and places it in out_path
        If out_path not given then just put it in current directory with the same basename

        @param readpaths - Path to the read to trim .fastq and .sff support only. Can be a tuple of paired reads
        @param qual_th - Quality threshold to trim reads on
        @param out_paths - Where to put the trimmed file[s]
        @param head_crop - How many bases to trim off the front

        @returns path to the trimmed fastq file
    '''
    readpaths = args[0]
    if isinstance(readpaths,str):
        readpaths = [readpaths]
    else:
        readpaths = sorted(readpaths)
    qual_th = args[1]
    if len(args) == 3:
        out_paths = args[2]
        if isinstance(out_paths,str):
            out_paths = [out_paths]
    else:
        out_paths = (None,None)
    headcrop = kwargs.get('head_crop', 0)
    primer_info = kwargs.get('primer_info')

    from Bio import SeqIO
    tfile = None
    for i, out_path in enumerate(out_paths):
        if out_path is None:
            out_paths[i] = basename( readpaths[i] ).replace('.sff','.fastq')
        logger.debug( "Using {0} as the output path".format(out_path) )

    # Keep the original name for later( Have to copy otherwise we are dealing with a pointer )
    orig_readpaths = [f for f in readpaths]

    # Convert sff to fastq
    for i,readpath in enumerate(readpaths):
        if readpath.endswith('.sff'):
            logger.debug( "Converting {0} to fastq".format(readpath) )
            # Just put in temp location then remove later
            _, tfile = tempfile.mkstemp(prefix='trimreads',suffix='sff.fastq')
            try:
                # Clip adapter based on clip_qual values in sff
                nwritten = reads.sffs_to_fastq( [readpath], tfile, True )
            except AssertionError as e:
                # Ignore the stupid sff bug in BioPython
                pass
            readpaths[i] = tfile

    # Run trimmer on the file
    trim_stats_dir = join( dirname(dirname(out_paths[0])), 'trim_stats' )
    stats_file = join( trim_stats_dir, basename(orig_readpaths[0]) + '.trim_stats' )
    if not isdir(dirname(stats_file)):
        os.makedirs( dirname(stats_file) )

    #run_cutadapt( readpath, stats=stats_file, o=out_path, q=qual_th )
    retpaths = []
    if len(readpaths) == 1:
        retpaths = [out_paths[0]]
        output = run_trimmomatic(
            'SE', readpaths[0], out_paths[0],
            ('LEADING',qual_th), ('TRAILING',qual_th), ('HEADCROP',headcrop),
            threads=1, trimlog=stats_file
        )
    else:
        retpaths = [out_paths[0],out_paths[0]+'.unpaired',out_paths[1],out_paths[1]+'.unpaired']
        output = run_trimmomatic(
            'PE', readpaths[0], readpaths[1], out_paths[0], out_paths[0]+'.unpaired', out_paths[1], out_paths[1]+'.unpaired',
            ('LEADING',qual_th), ('TRAILING',qual_th), ('HEADCROP',headcrop),
            threads=1, trimlog=stats_file, primer_info=primer_info #NOTE: only run primer on paired read files
        )

    # Prepend stats file with stdout from trimmomatic
    with open(stats_file, 'r+') as fh:
        contents = fh.read()
        fh.seek(0)
        fh.write( output )
        fh.write( '\n' )
        fh.write( contents )

    # Clean up temp file
    if tfile:
        os.unlink(tfile)

    return retpaths

def run_trimmomatic( *args, **kwargs ):
    '''
        Runs trimmomatic
        @param arg0 - SE|PE -- Only SE supported at this time
        @param arg1-arg6 - Input/Ouput files
        @param arg7-argN - Tuples of (Trimmomatic Step,Options)
        @param kwargs are any --options

        run_trimmomatic( 'SE', 'input.fq', 'output.fq', ('LEADING','20), ('TRAILING','20'), trim_log='out.log' )
        would result in
        java -jar trimmomatic.jar input.fq output.fq LEADING:20 TRAILING:20 --trim_log out.log
    '''
    if args[0] == 'SE':
        inputs = [args[1]]
        outputs = [args[2]]
        # Trimmomatic doesn't seem to be able to detect Sanger quality encoding
        # so we will try to force it here to phred33
        if data.is_sanger_readfile(args[1]):
            kwargs['phred33'] = ''
        steps = args[3:]
    elif args[0] == 'PE':
        inputs = list(args[1:3])
        outputs = list(args[3:7])
        steps = args[7:]
    else:
        raise ValueError( 'SE or PE need to be supplied' )

    # Change all steps to strings of STEPNAME:VALUE
    steps = [':'.join([str(x) for x in s]) for s in steps]
    # Set all options
    options = shlex.split( ' '.join( ['-{0} {1}'.format(k,v) for k,v in kwargs.items()] ) )
    cmd = ['trimmomatic', args[0]] + options + inputs + outputs + steps

    primer_info = kwargs.get('primer_info')
    if primer_info:
        cmd += ':'.join(['ILLUMINACLIP'] + primer_info)


    # Write stdout to output argument(should be fastq)
    # Allow us to read stderr which should be stats from cutadapt
    logger.debug( "Running {0}".format(' '.join(cmd)) )
    try:
        output = compat.check_output( cmd, stderr=subprocess.STDOUT )
        return output
    except subprocess.CalledProcessError as e:
        logger.critical( "Trimmomatic error: {0}".format(e.output) )
        raise e

def run_cutadapt( *args, **kwargs ):
    '''
        Runs cutadapt with the given arguments and kwargs

        @param - fastq file to trim
        @param - output file location
        @param q - Quality threshold

        @returns the stderr output from cutadapt
    '''
    outpath = kwargs.get('o')
    cmd = ['cutadapt', '-o', outpath, '-q', str(kwargs.get('q')), args[0]]
    out_stats = kwargs.get( 'stats', outpath + '.trim_stats' )
    fout = open(out_stats,'wb')
    # Write stdout to output argument(should be fastq)
    # Allow us to read stderr which should be stats from cutadapt
    logger.debug( "Running {0}".format(cmd) )
    logger.debug( "Sending stdout to {0}".format(out_stats) )
    p = subprocess.Popen( cmd, stdout=fout, stderr=subprocess.PIPE )
    # Only stderr should be available
    _,se = p.communicate()
    if p.returncode != 0:
        e = subprocess.CalledProcessError(p.returncode,' '.join(cmd), se)
        raise e
    return se

def parse_args( args=sys.argv[1:] ):
    from ngs_mapper import config
    conf_parser, args, config, configfile = config.get_config_argparse(args)
    defaults = config['trim_reads']

    parser = argparse.ArgumentParser(
        parents=[conf_parser],
        description='Trims reads'
    )

    parser.add_argument(
        dest='readsdir',
        help='Read or directory of read files'
    )

    parser.add_argument(
        '-q',
        dest='q',
        default=defaults['q']['default'],
        help=defaults['q']['help']
    )

    parser.add_argument(
        '--head-crop',
        dest='headcrop',
        default=defaults['headcrop']['default'],
        help=defaults['headcrop']['help']
    )

    parser.add_argument(
        '-o',
        dest='outputdir',
        default=defaults['outputdir']['default'],
        help=defaults['outputdir']['help']
    )

    parser.add_argument(
        '--platforms',
        dest='platforms',
        choices=defaults['platforms']['choices'],
        default=defaults['platforms']['default'],
        help=defaults['platforms']['help']
    )

    parser.add_argument(
        '--primer-file',
        dest='primer_file',
        default=defaults['primerfile']['default'],
        help=defaults['primerfile']['help']
    )

    parser.add_argument(
        '--primer-seed',
        dest='primer_seed',
        default=defaults['primerseed']['default'],
        help=defaults['primerseed']['help']
    )

    parser.add_argument(
        '--palindrome-clip',
        dest='palindrom_clip',
        default=defaults['palindromeclip']['default'],
        help=defaults['palindromeclip']['help']
    )

    parser.add_argument(
        '--simple-clip',
        dest='simple_clip',
        default=defaults['simpleclip']['default'],
        help=defaults['simpleclip']['help']
    )

    return parser.parse_args( args )
