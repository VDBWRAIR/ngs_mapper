"""
The intention of this script is to easily sync a given IonTorrent or IonProton run path into the :doc:`NGSData <../ngsdata>` structure

You need to ensure that the run directory for the run you want to sync is available and that you know the path to it.

The run directores for IonTorrent and IonProton are located under /data/analysis/Home on the instrument.
These run directories should contain a ion_params_00.json as well as plugin_out/downloads that contains the fastq files for each barcode.

Usage
=====

.. code-block:: bash

    ion_sync /path/to/run_directory

Optionally you can specify the NGSdata path as follows

.. code-block:: bash

    ion_sync /path/to/run_directory --ngsdata /path/to/NGSData

How it works
============

#. The run_directory is copied into the RawData directory located under the specified ``--ngsdata``
#. Then a sample mapping is created based on the ``ion_params_00.json`` file that is located inside of the run directory.

    * The ion_params_00.json is automatically created by the instrument based 
      on the information provided in the web interface for the run.
    * The ``ion_params_00.json`` file is loaded and the samplenames are taken from
      experimentAnalysisSettings->barcodedSamples inside of the loaded file.

#. All of the fastq files inside of plugin_out/downloads are compiled and all 
   IonXpress_XXX are symlinked into the ReadData under the path specified with
   ``--ngsdata`` and renamed using the mapping that was created from 
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

Additionally, the naming format of the fastq files as well as the placement of them may change with the different versions of the software.
If you notice that your fastq files are not located inside of ``plugin_out/downloads`` inside of your runs, please :doc:`submit an issue <../createissue>` and provide the relative path to your run that your fastq files are located.

If you only have .bam files you should be able to convert them to fastq with the following command

.. code-block:: bash

    samtools view /path/to/IonXpress_XXX.bam | awk '{printf("@%s\\n%s\\n+\\n%s\\n",$1,$10,$11)}' > /path/to/IonXpress_XXX.fastq

"""

import json
from os.path import dirname, basename, join, exists, isdir, relpath
import argparse
import glob
import shutil
import sys
import os
import argparse

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

def get_samplefile_mapping(barcodemapping, fastqs):
    '''
    Given a barcode mapping and a list of fastq paths, returns the mapping
    that can be used to rename those files
    The returned mapping will be keyed by the original file path given and 
    the key will be just the new filename without the path
    
    If the barcode is missing from the barcodemapping then just 
    the filename will be used for the value

    :param dict barcodemapping: {'barcode1':'samplename1',...}
    :param list fastqs: list of fastq barcode named files. Probably plugin_out/downloads/IonXpress_001.whatever.fastq
    :return: mapping of {'originalfilepath:'newfilename no path'}
    '''
    rename_map = {}
    for fq in fastqs:
        path = dirname(fq)
        filename = basename(fq)
        p = filename.split('.')
        if len(p) == 3:
            barcode, run, ext = p
        else:
            raise InvalidFastqFilename("{0} is an invalid filename. Contains more than 2 '.' in the filename".format(fq))
        try:
            newname = barcodemapping[barcode]
        except KeyError as e:
            newname = barcode
        rename_map[fq] = '{0}.{1}.{2}.{3}'.format(newname,barcode,run,ext)
    return rename_map

def ion_mapping(fastqpath, ionparampath):
    '''
    Get a mapping of original filepath to new filename given
    the path to the fastqs and ion_param_00.json file

    :param str fastqpath: Path to fastq directory
    :param str ionparmpath: path to ion_param*.json
    :return: dictionary of {'original path': 'newpath with samplename'}
    '''
    fastqs = glob.glob(join(fastqpath,'*.fastq'))
    ionparam = json.load(open(ionparampath))
    samplemap = get_samplemapping(ionparam)
    samplefilemap = get_samplefile_mapping(samplemap, fastqs)
    return samplefilemap

def sync_run(runpath, ngsdata, printmappingonly):
    '''
    Sync an iontorrent or ionproton run into ngsdata
    The runs need to contain ion_params_00.json as well as
    plugin_out/downloads/IonXpress_XXX.runname.fastq

    :param str runpath: path to ion run
    :param str ngsdata: path to ngsdata
    :param bool printmappingonly: True to only print mapping
    '''
    runname = basename(runpath)
    rawdata = join(ngsdata, 'RawData', 'IonTorrent', runname)
    readdata = join(ngsdata, 'ReadData', 'IonTorrent', runname)
    rbs = join(ngsdata, 'ReadsBySample')
    # Have to give relative path to fastq files later on
    cwd = os.getcwd()
    os.chdir(runpath)
    fastqpath = join('plugin_out','downloads')
    ionparampath = 'ion_params_00.json'
    # This could raise an exception
    samplefilemap = ion_mapping(fastqpath, ionparampath)
    # Ensure we return to the cwd
    os.chdir(cwd)

    if printmappingonly:
        for k,v in samplefilemap.items():
            sys.stdout.write("{0} -> {1}\n".format(k,v))
        return

    if isdir(rawdata):
        sys.stderr.write(
            '{0} already exists so will not be synced.\n' \
            'You may need to remove this directory and rerun the '\
            'sync\n'.format(rawdata)
        )
    else:
        try:
            shutil.copytree(runpath, rawdata)
        except shutil.Error as e:
            # Probably doesn't matter because broken symlink
            pass
    sync_readdata(samplefilemap, readdata)
    sync_readsbysample(readdata, rbs)

def sync_readdata(samplefilemap, readdatapath):
    '''
    Sync a list of raw fastq files into the ReadData

    :param dict samplefilemap: mapping of origin path to newname
    :param str readdatapath: path to sync into
    '''
    runname = basename(readdatapath)
    sympathbase = join('../../../RawData/IonTorrent', runname)
    if not isdir(readdatapath):
        os.makedirs(readdatapath)
    for origpath, newname in samplefilemap.iteritems():
        src = join(sympathbase, origpath)
        dst = join(readdatapath, newname)
        if exists(dst):
            sys.stderr.write('{0} exists and will be skipped\n'.format(dst))
            continue
        os.symlink(src, dst)

def sync_readsbysample(readdatapath, readsbysample):
    '''
    Sync the ReadData into the readsbysample

    :param str readdatapath: path to readsbysample for run being synced
    :param str readsbysample: path to readsbysample
    '''
    # Relative path between ReadsBysample and ReadData(../)
    sympathbase = relpath(readdatapath, readsbysample)
    for fq in os.listdir(readdatapath):
        parts = fq.split('.')
        # Have to add another level deeper for relative symlink
        src = join('..', sympathbase, fq)
        sn = parts[0]
        snpath = join(readsbysample, sn)
        if not isdir(snpath):
            os.makedirs(snpath)
        dst = join(snpath, fq)
        if exists(dst):
            sys.stderr.write('{0} already exists\n'.format(dst))
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

    return parser.parse_args(args)

def main():
    args = parse_args()
    sync_run(args.rundir, args.ngsdata, args.print_samplemapping)
