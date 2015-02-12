"""
The intention of this script is to easily sync a given IonTorrent or IonProton run path into the :doc:`NGSData <../ngsdata>` structure

You need to ensure that the run directory for the run you want to sync is available and that you know the path to it.

The run directores for IonTorrent and IonProton are located under /data/analysis/Home on the instrument.
These run directories should contain a ion_params_00.json as well as plugin_out/downloads that contains the fastq files for each barcode.

Usage
=====

.. code-block:: bash

    $> ion_sync /path/to/run_directory

Optionally you can specify the NGSdata path as follows

.. code-block:: bash

    $> ion_sync /path/to/run_directory --ngsdata /path/to/NGSData

Preview the samplemapping
=========================

You can see how the IonXpress_XXX barcoded files will be renamed by using the
``--pring-samplemapping`` option

.. code-block:: bash

    $> ion_sync --print-samplemapping run_1/
    run_1/basecaller_results/IonXpress_002_rawlib.basecaller.bam -> sample2.IonXpress_002.run_1.fastq
    run_1/basecaller_results/IonXpress_003_rawlib.basecaller.bam -> IonXpress_003.IonXpress_003.run_1.fastq
    run_1/basecaller_results/IonXpress_001_rawlib.basecaller.bam -> sample1.IonXpress_001.run_1.fastq

How it works
============

#. The run_directory is copied into the RawData directory located under the specified ``--ngsdata``
#. Then a sample mapping is created based on the ``ion_params_00.json`` file that is located inside of the run directory.

    * The ion_params_00.json is automatically created by the instrument based 
      on the information provided in the web interface for the run.
    * The ``ion_params_00.json`` file is loaded and the samplenames are taken from
      experimentAnalysisSettings->barcodedSamples inside of the loaded file.

#. All of the fastq files inside of plugin_out/downloads are compiled and all 
#. All IonXpress_XXX are symlinked into the ReadData under the path specified with
   ``--ngsdata`` if their filesize is greater than the --min-fastq-size
#. Then they are renamed using the mapping that was created from 
   the ``ion_params_00.json`` 
#. All of the fastq files under ReadData are then symlinked under their 
   cooresponding samplename directory inside of ReadsBySample under the path
   specified by ``--ngsdata``

    * samplenames are determined by splitting the fastq file name by each 
      occurance of a dot and taking the first item from that.

Notes
=====

As the IonTorrent's and IonProton's software versions change sometimes so does the way the ``ion_params_00.json`` is created by the instrument.
This means that you may encounter an error as the location of the samplename -> barcode mapping inside of the ``ion_params_00.json`` may
change. As with any issue, if you encounter an error while using this script, :doc:`submit an issue <../createissue>` and make sure to copy-paste the contents of your ion_params_00.json file in the issue.
"""

import json
import os.path
import argparse
import glob
import shutil
import sys
import os
import argparse
import re
import bam
import log

logger = log.setup_logger(__name__, log.get_config())

# MegaByte size-ish
MB = 2**20

from ngs_mapper import config

class InvalidIonParam(Exception): pass

def get_samplemapping(ionparam):
    '''
    :param json ionparam: json parsed ionparam
    :return: sample mapping of {'barcode': 'samplename', ...}
    '''
    try:
        samples_str = ionparam['experimentAnalysisSettings']['barcodedSamples']
    except KeyError as e:
        raise InvalidIonParam(str(e))

    if not samples_str:
        samples_str = '{}'
    samples_json = json.loads(samples_str)
    mapping = {'nomatch': 'nomatch'}
    for sample, info in samples_json.iteritems():
        for barcode in info['barcodes']:
            mapping[barcode] = sample
    return mapping

class InvalidFastqFilename(Exception): pass

def get_samplefile_mapping(barcodemapping, fastqs, runname):
    '''
    Given a barcode mapping and a list of fastq|rawlib.basecaller.bam paths, 
    returns the mapping that can be used to rename those files
    The returned mapping will be keyed by the original file path given and 
    the key will be just the new filename without the path

    Possible examples:
    path/to/IonXpress_XXX_rawlib.basecaller.bam
    path/to/IonXpress_XXX.run_name.fastq
    
    If the barcode is missing from the barcodemapping then just 
    the filename will be used for the value

    :param dict barcodemapping: {'barcode1':'samplename1',...}
    :param list fastqs: list of fastq|bam barcode named files
    :param str runname: what to put in the 3rd field of fastq name to make unique
    :return: mapping of {'originalfilepath:'newfilename no path'}
    '''
    rename_map = {}
    for fq in fastqs:
        path = os.path.dirname(fq)
        filename = os.path.basename(fq)
        m = re.search('(IonXpress_\d+|nomatch)(?:_rawlib){0,1}.(\w+).(fastq|bam)', filename)
        if m:
            barcode, run, ext = m.groups()
        else:
            raise InvalidFastqFilename("{0} is an invalid filename".format(fq))
        try:
            newname = barcodemapping[barcode]
        except KeyError as e:
            newname = barcode
        rename_map[fq] = '{0}.{1}.{2}.{3}'.format(newname,barcode,runname,'fastq')
    return rename_map

def ion_mapping(reads, ionparampath):
    '''
    Get a mapping of original filepath to new filename given
    the path to the fastqs and ion_param_00.json file

    :param str reads: Path to directory of fastq or bam files
    :param str ionparmpath: path to ion_param*.json
    :return: dictionary of {'original path': 'newpath with samplename'}
    '''
    #reads = glob.glob(os.path.join(fastqpath,'{IonXpress,nomatch}*{bam,fastq}'))
    ionparam = json.load(open(ionparampath))
    expName = ionparam.get('expName', '')
    if not expName:
        raise InvalidIonParam('{0} is missing expName'.format(ionparampath))
    runname = ionparam['expName']
    samplemap = get_samplemapping(ionparam)
    samplefilemap = get_samplefile_mapping(samplemap, reads, runname)
    return samplefilemap

def convert_basecaller_results_to_fastq(bamfastqmapping, fastqoutpath):
    '''
    Convert list of bam files into fastq using the supplied mapping

    Ensures that fastqoutpath exists

    :param dict bamfastqmapping: mapping of rawlibbams to fastq basenames
    :param str fastqoutpath: directory to place converted fastq in
    '''
    if not os.path.exists(fastqoutpath):
        logger.debug('Creating {0} to put output from bam conversion in'.format(
            fastqoutpath
        ))
        os.makedirs(fastqoutpath)
    for bamf, fqf in bamfastqmapping.iteritems():
        outfq_path = os.path.join(fastqoutpath, fqf)
        with open(outfq_path, 'w') as fh:
            logger.debug('Converting {0} to {1}'.format(bamf, outfq_path))
            for read in bam.bam_to_fastq(bamf):
                fh.write(read + '\n')

def sync_run(runpath, ngsdata, printmappingonly, minfastqsize):
    '''
    Sync an iontorrent or ionproton run into ngsdata
    The runs need to contain ion_params_00.json as well as
    plugin_out/downloads/IonXpress_XXX.runname.fastq

    :param str runpath: path to ion run
    :param str ngsdata: path to ngsdata
    :param bool printmappingonly: True to only print mapping
    '''
    runname = os.path.basename(runpath)
    rawdata = os.path.join(ngsdata, 'RawData', 'IonTorrent', runname)
    readdata = os.path.join(ngsdata, 'ReadData', 'IonTorrent', runname)
    rbs = os.path.join(ngsdata, 'ReadsBySample')
    logger.debug("Run Name: {0}\nRawData: {1}\nReadData: {2}\nReadsBySample: {3}\n".format(
        runname, rawdata, readdata, rbs
    ))
    # Have to give relative path to fastq files later on
    #cwd = os.getcwd()
    #os.chdir(runpath)
    ionparampath = os.path.join(runpath, 'ion_params_00.json')
    basecaller_results_dir = os.path.join(runpath, 'basecaller_results')
    logger.debug('Looking for rawlib bam files inside of {0}'.format(
        basecaller_results_dir
    ))
    rawlibbams = glob.glob(os.path.join(basecaller_results_dir,'*_rawlib*.bam'))
    logger.debug('Found the following rawlib bam files: {0}'.format(rawlibbams))
    bamfastqmap = ion_mapping(rawlibbams, ionparampath)
    # Ensure we return to the cwd
    #os.chdir(cwd)

    if printmappingonly:
        for k,v in bamfastqmap.items():
            sys.stdout.write("{0} -> {1}\n".format(k,v))
        return

    if os.path.isdir(rawdata):
        logger.info(
            '{0} already exists so the raw data directory will not be copied again.\n' \
            'You may need to remove this directory and rerun the '\
            'sync\n'.format(rawdata)
        )
    else:
        try:
            shutil.copytree(runpath, rawdata)
        except shutil.Error as e:
            # Probably doesn't matter because broken symlink
            pass

    # reset to the copied rawdata version
    ionparampath = os.path.join(rawdata, 'ion_params_00.json')
    fastqpath = os.path.join(rawdata, 'plugin_out', 'downloads')
    # If no plugin_out/downloads, make it from basecaller_results/bams
    if not glob.glob(os.path.join(fastqpath,'IonXpress_*.fastq')):
        bamfqmap = {}
        logger.info('Putting converted bam->fastq into {0}'.format(fastqpath))
        logger.debug('Basecaller_results bam -> fastq mapping')
        for b, f in bamfastqmap.iteritems():
            # Strip off the sample name off front
            bamfqmap[b] = '.'.join(f.split('.')[1:])
            logger.debug('{0} -> {1}'.format(b,bamfqmap[b]))
        convert_basecaller_results_to_fastq(bamfqmap, fastqpath)
    # Fastqs should be guaranteed to exist now...hopefully
    fastqs = [os.path.join('plugin_out', 'downloads',fq) for fq in os.listdir(fastqpath)]
    logger.debug('Found the following fastq files: {0}'.format(fastqs))
    # This could raise an exception
    samplefilemap = ion_mapping(fastqs, ionparampath)
    logger.debug('Sample fastq file mapping:')
    for o, n in samplefilemap.iteritems():
        logger.debug('{0} -> {1}'.format(o,n))
    logger.info('Syncing read data')
    sync_readdata(samplefilemap, readdata, minfastqsize)
    logger.info('Syncing reads by sample')
    sync_readsbysample(readdata, rbs)

def sync_readdata(samplefilemap, readdatapath, minfastqsize):
    '''
    Sync a list of raw fastq files into the ReadData

    :param dict samplefilemap: mapping of origin path to newname
    :param str readdatapath: path to sync into
    '''
    runname = os.path.basename(readdatapath)
    sympathbase = os.path.join('../../../RawData/IonTorrent', runname)
    if not os.path.isdir(readdatapath):
        logger.debug('Creating {0} to put fastq files'.format(readdatapath))
        os.makedirs(readdatapath)
    for origpath, newname in samplefilemap.iteritems():
        src = os.path.join(sympathbase, origpath)
        dst = os.path.join(readdatapath, newname)
        _origpath = os.path.join(readdatapath, src)
        if os.path.exists(dst):
            logger.info('{0} exists and will be skipped'.format(dst))
            continue
        filesize = os.stat(_origpath).st_size
        if filesize < minfastqsize:
            logger.info(
                '{0} will be skipped as its filesize({1}) is below the min-fastq-size ' \
                'set({2})'.format(
                    _origpath, filesize, minfastqsize
            ))
            continue
        logger.debug('Symlinking {0} -> {1}'.format(src,dst))
        os.symlink(src, dst)

def sync_readsbysample(readdatapath, readsbysample):
    '''
    Sync the ReadData into the readsbysample

    :param str readdatapath: path to readsbysample for run being synced
    :param str readsbysample: path to readsbysample
    '''
    # Relative path between ReadsBysample and ReadData(../)
    sympathbase = os.path.relpath(readdatapath, readsbysample)
    for fq in os.listdir(readdatapath):
        parts = fq.split('.')
        # Have to add another level deeper for relative symlink
        src = os.path.join('..', sympathbase, fq)
        sn = parts[0]
        snpath = os.path.join(readsbysample, sn)
        if not os.path.isdir(snpath):
            os.makedirs(snpath)
        dst = os.path.join(snpath, fq)
        if os.path.exists(dst):
            logger.info('{0} already exists'.format(dst))
            continue
        os.symlink(src, dst)

def parse_args():
    conf_parser, args, cfg, configfile = config.get_config_argparse(sys.argv[1:])
    defaults = cfg['ion_sync']

    parser = argparse.ArgumentParser(
        description='Sync roche run directories into the NGSData data structure',
        parents=[conf_parser]
    )

    parser.add_argument(
        'rundir',
        help='Path to IonTorrent or IonProton run which contains' \
            ' an ion_params_00.json file and a ' \
            ' plugin_out/downloads directory with fastq files'
    )

    parser.add_argument(
        '--ngsdata',
        default=defaults['ngsdata']['default'],
        help=defaults['ngsdata']['help']
    )

    parser.add_argument(
        '--print-samplemapping',
        dest='print_samplemapping',
        default=False,
        action='store_true',
        help='Just print the sample mapping for the run and quit'
    )

    parser.add_argument(
        '--min-fastq-size',
        dest='min_fastq_size',
        type=int,
        default=defaults['min_fastq_size']['default'],
        help=defaults['min_fastq_size']['help']
    )

    return parser.parse_args(args)

def main():
    args = parse_args()
    sync_run(args.rundir, args.ngsdata, args.print_samplemapping, args.min_fastq_size)
