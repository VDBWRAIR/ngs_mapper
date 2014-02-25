#!/usr/bin/env python

import argparse
import subprocess
import shlex
import sys
import os
import os.path
import tempfile
import logging
import shutil
import glob

# Everything to do with running a single sample
# Geared towards running in a Grid like universe(HTCondor...)
# Ideally the entire sample would be run inside of a prefix directory under
# /dev/shm and drop back on tmpdir if /dev/shm didn't exist

logconfig= logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)-15s %(message)s'
)
log = logging.getLogger(logconfig)

class MissingCommand(Exception):
    pass

class AlreadyExists(Exception):
    pass

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description='Runs a single sample through the pipeline'
    )

    parser.add_argument(
        dest='readsdir',
        help='Directory that contains reads to be mapped'
    )

    parser.add_argument(
        dest='reference',
        help='The path to the reference to map to'
    )

    parser.add_argument(
        dest='prefix',
        help='The prefix to put before every output file generated. Probably the samplename'
    )

    default_outdir = os.getcwd()
    parser.add_argument(
        '-od',
        '--outdir',
        dest='outdir',
        default=default_outdir,
        help='The output directory for all files to be put[Default: {}]'.format(default_outdir)
    )

    return parser.parse_args( args )

def temp_projdir( prefix, suffix='runsample' ):
    '''
        Get a temporary directory that all files for the sample can be compiled into
        Prefer /dev/shm

        @param prefix/suffix - look these up in tempdir.mkdtemp docs

        @returns the temporary directory used
    '''
    tdir = '/dev/shm'
    if not os.path.isdir( '/dev/shm' ):
        tdir = '/tmp'

    return tempfile.mkdtemp( suffix, prefix, dir=tdir )

def run_cmd( cmdstr, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr, script_dir=os.path.dirname(__file__) ):
    '''
        Runs a subprocess on cmdstr and logs some timing information for each command

        @params stdin/out/err should be whatever is acceptable to subprocess.Popen for the same
        
        @returns the popen object
    '''
    cmd = shlex.split( cmdstr )
    cmd[0] = os.path.join( script_dir, cmd[0] )
    log.debug( "Running {}".format(' '.join(cmd)) )
    try:
        p = subprocess.Popen( cmd, stdout=stdout, stderr=stderr, stdin=stdin )
        return p
    except OSError as e:
        raise MissingCommand( "{} is not an executable?".format(cmd[0]) )

def main( args ):
    tdir = temp_projdir( args.prefix )
    bamfile = os.path.join( tdir, args.prefix + '.bam' )
    flagstats = os.path.join( tdir, 'flagstats.txt' )
    consensus = os.path.join( tdir, bamfile+'.consensus.fastq' )
    bwalog = os.path.join( tdir, 'bwa.log' )
    logfile = os.path.join( tdir, args.prefix + '.log' )

    if os.path.isdir( args.outdir ):
        if os.listdir( args.outdir ):
            raise AlreadyExists( "{} already exists and is not empty".format(args.outdir) )

    log.debug( "--- Starting {} --- ".format(args.prefix) )
    with open(logfile,'wb') as lfile:
        cmd_args = {
            'tdir': tdir,
            'readsdir': args.readsdir,
            'reference': args.reference,
            'bamfile': bamfile,
            'flagstats': flagstats,
            'consensus': consensus
        }

        # Best not to run across multiple cpu/core/threads on any of the pipeline steps
        # as multiple samples may be running concurrently already

        # Return code list
        rets = []

        # Mapping
        with open(bwalog, 'wb') as blog:
            cmd = 'run_bwa_on_samplename.py {readsdir} {reference} -o {bamfile} -t 1'
            p = run_cmd( cmd.format(**cmd_args), stdout=blog, stderr=subprocess.STDOUT )
            # Wait for the sample to map
            rets.append( p.wait() )
            # Everything else is dependant on bwa finishing so might as well die here
            if rets[-1] != 0:
                cmd = cmd.format(**cmd_args)
                log.critical( "{} failed to complete sucessfully. Please check the log file {} for more details".format(cmd,bwalog) )
                sys.exit(1)

        # Tag Reads
        cmd = 'tagreads.py {bamfile}'
        p = run_cmd( cmd.format(**cmd_args), stdout=lfile, stderr=subprocess.STDOUT )
        rets.append( p.wait() )

        # Variant Calling
        cmd = 'base_caller.py {bamfile} {reference}'
        p = run_cmd( cmd.format(**cmd_args), stdout=lfile, stderr=subprocess.STDOUT )
        rets.append( p.wait() )
        if rets[-1] != 0:
            cmd = cmd.format(**cmd_args)
            log.critical( '{} failed to complete successfully'.format(cmd) )

        # Flagstats
        with open(flagstats,'wb') as flagstats:
            cmd = 'samtools flagstat {bamfile}'
            p = run_cmd( cmd.format(**cmd_args), stdout=flagstats, stderr=lfile, script_dir='' )
            rets.append( p.wait() )

        # Graphics
        cmd = 'graphsample.py {bamfile} -od {tdir}'
        p = run_cmd( cmd.format(**cmd_args), stdout=lfile, stderr=subprocess.STDOUT )
        rets.append( p.wait() )

        # Consensus
        with open(consensus,'wb') as consensus:
            cmd = 'gen_consensus.sh {reference} {bamfile}'
            p = run_cmd( cmd.format(**cmd_args), stdout=consensus, stderr=lfile )
            rets.append( p.wait() )

        # If sum is > 0 then one of the commands failed
        if sum(rets) != 0:
            log.critical( "!!! There was an error running part of the pipeline !!!" )
            log.critical( "Please check the logfile {}".format(logfile) )
            sys.exit( 1 )

        log.debug( "Moving {} to {}".format( tdir, args.outdir ) )
        if not os.path.isdir( args.outdir ):
            shutil.move( tdir, args.outdir )
        else:
            file_list = glob.glob( os.path.join( tdir, '*' ) )
            for f in file_list:
                shutil.move( f, args.outdir )
        log.debug( "--- Finished {} ---".format(args.prefix) )

if __name__ == '__main__':
    main(parse_args())
