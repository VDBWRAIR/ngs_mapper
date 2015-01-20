#!/usr/bin/env python

import json
from os.path import dirname, basename, join
import argparse
import glob

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
