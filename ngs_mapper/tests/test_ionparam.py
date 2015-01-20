import json

from unittest2 import TestCase
import mock

from .. import ionparam

class Base(TestCase):
    def setUp(self):
        self.barcoded_samples = {
            'sample1': {
                'barcodes': ['IonXpress_001']
            },
            'sample2': {
                'barcodes': ['IonXpress_002', 'IonXpress_003']
            }
        }

        self.ionparam = {
            'experimentAnalysisSettings': {
                'barcodedSamples': json.dumps(self.barcoded_samples)
            }
        }

        # _004 purposely has no matching mapping in ionparam
        self.fastqs = [
            'plugin_out/downloads/IonXpress_001.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_002.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_003.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_004.R_run_name.fastq',
            'plugin_out/downloads/nomatch.R_run_name.fastq',
        ]

        self.samplemap = {
            'IonXpress_001': 'sample1',
            'IonXpress_002': 'sample2',
            'IonXpress_003': 'sample2',
            # Missing IonXpress_004
            'nomatch': 'nomatch',
        }

        self.samplefilemap = {
            'plugin_out/downloads/IonXpress_001.R_run_name.fastq':
                'sample1.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_002.R_run_name.fastq':                 'sample2.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_003.R_run_name.fastq':
                'sample2.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_004.R_run_name.fastq':
                'IonXpress_004.R_run_name.fastq',
            'plugin_out/downloads/nomatch.R_run_name.fastq':
                'nomatch.R_run_name.fastq',
        }

class TestGetSamplemapping(Base):
    def test_ionparam_missing_keys(self):
        self.assertRaises(
            ionparam.InvalidIonParam,
            ionparam.get_samplemapping, {}
        )

    def test_empty_barcoded_samples(self):
        self.ionparam['experimentAnalysisSettings']['barcodedSamples'] = ''
        self.assertEqual(
            {'nomatch': 'nomatch'},
            ionparam.get_samplemapping(self.ionparam)
        )

    def test_creates_correct_mapping(self):
        expect = self.samplemap
        r = ionparam.get_samplemapping(self.ionparam)
        for bc in expect:
            self.assertEqual(
                expect[bc], r[bc]
            )

class TestGetSamplefileMapping(Base):
    def test_creates_correct_mapping(self):
        expect = self.samplefilemap
        r = ionparam.get_samplefile_mapping(
            self.samplemap,
            self.fastqs
        )
        for e in expect:
            self.assertEqual(
                expect[e],
                r[e]
            )

    def test_fastq_name_more_than_2_periods(self):
        self.fastqs = ['IonXpress_001.R.run_name.fastq']
        self.assertRaises(
            ionparam.InvalidFastqFilename,
            ionparam.get_samplefile_mapping,
            self.samplemap,
            self.fastqs
        )

@mock.patch.object(ionparam, 'glob')
@mock.patch.object(ionparam, 'json')
@mock.patch('__builtin__.open', mock.Mock())
class TestIonMapping(Base):
    def test_creates_correct_mapping(self, mjson, mglob):
        mjson.load.return_value = self.ionparam
        mjson.loads.return_value = self.barcoded_samples
        mglob.glob.return_value = self.fastqs
        r = ionparam.ion_mapping(
            'plugin_out/downloads',
            'ion_param_00.json'
        )
        for origpath in self.samplefilemap:
            self.assertEqual(
                self.samplefilemap[origpath],
                r[origpath]
            )
