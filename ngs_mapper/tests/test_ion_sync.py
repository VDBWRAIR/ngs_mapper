import tempfile
import os
from os.path import *
import shutil
from glob import glob
import re
import json
from StringIO import StringIO

import mock
from unittest2 import TestCase

from .. import ion_sync

FILE = mock.Mock()

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
                'sample1.IonXpress_001.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_002.R_run_name.fastq':                 'sample2.IonXpress_002.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_003.R_run_name.fastq':
                'sample2.IonXpress_003.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_004.R_run_name.fastq':
                'IonXpress_004.IonXpress_004.R_run_name.fastq',
            'plugin_out/downloads/nomatch.R_run_name.fastq':
                'nomatch.nomatch.R_run_name.fastq',
        }

class TestGetSamplemapping(Base):
    def test_ionparam_missing_keys(self):
        self.assertRaises(
            ion_sync.InvalidIonParam,
            ion_sync.get_samplemapping, {}
        )

    def test_empty_barcoded_samples(self):
        self.ionparam['experimentAnalysisSettings']['barcodedSamples'] = ''
        self.assertEqual(
            {'nomatch': 'nomatch'},
            ion_sync.get_samplemapping(self.ionparam)
        )

    def test_creates_correct_mapping(self):
        expect = self.samplemap
        r = ion_sync.get_samplemapping(self.ionparam)
        for bc in expect:
            self.assertEqual(
                expect[bc], r[bc]
            )

class TestGetSamplefileMapping(Base):
    def test_creates_correct_mapping(self):
        expect = self.samplefilemap
        r = ion_sync.get_samplefile_mapping(
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
            ion_sync.InvalidFastqFilename,
            ion_sync.get_samplefile_mapping,
            self.samplemap,
            self.fastqs
        )

@mock.patch.object(ion_sync, 'glob')
@mock.patch.object(ion_sync, 'json')
@mock.patch('__builtin__.open', mock.Mock())
class TestIonMapping(Base):
    def test_creates_correct_mapping(self, mjson, mglob):
        mjson.load.return_value = self.ionparam
        mjson.loads.return_value = self.barcoded_samples
        mglob.glob.return_value = self.fastqs
        r = ion_sync.ion_mapping(
            'plugin_out/downloads',
            'ion_param_00.json'
        )
        for origpath in self.samplefilemap:
            self.assertEqual(
                self.samplefilemap[origpath],
                r[origpath]
            )

class BaseSync(Base):
    def setUp(self):
        super(BaseSync,self).setUp()
        self.tdir = tempfile.mkdtemp(prefix='ionsync',suffix='tests')
        self.ngsdata = join(self.tdir,'NGSDATA')
        self.rawdata = join(self.ngsdata,'RawData')
        self.readdata = join(self.ngsdata,'ReadData')
        self.readsbysample = self.rbs = join(self.ngsdata,'ReadsBySample')
        os.chdir(self.tdir)

    def tearDown(self):
        os.chdir('/')
        shutil.rmtree(self.tdir)

    def _create_rundir(self, reads, ionparam, path=None):
        '''
        Create ion[Torrent|Proton] run directory
        reads is a list of IonXpress_XXX.whatever.fastq
        ionparam is a dict representing the json.load from actual
        ion_params_00.json file
        '''
        if path is None:
            rundir = tempfile.mkdtemp(prefix='rundir_',suffix='ion',dir=self.tdir)
        else:
            rundir = path
        downloads = join(rundir,'plugin_out','downloads')
        ion_params = join(rundir,'ion_params_00.json')
        # Make plugin_out/downloads dir
        os.makedirs(downloads)
        # Make ion_params_00.json file
        with open(ion_params,'w') as fh:
            json.dump(ionparam, fh)
        # Make mock fastqs
        for read in reads:
            readpath = join(rundir,read)
            open(readpath,'w').close()
        return rundir, downloads, ion_params

@mock.patch.object(ion_sync, 'config')
@mock.patch.object(ion_sync, 'argparse')
class TestFunctional(BaseSync):
    def setUp(self):
        super(TestFunctional,self).setUp()
        self.rundir, self.downloads, self.ion_params = \
            self._create_rundir(self.fastqs, self.ionparam)
        # Make a broken symlink to make sure
        # it is handled with shutil.copytree
        os.symlink(
            '/missing/path/IonXpress_001.bam',
            join(self.downloads,'IonXpress_001.bam')
        )
        self.runname = basename(self.rundir)
        self.readdata = join(self.readdata,'IonTorrent',self.runname)
        self.rawdata = join(self.rawdata,'IonTorrent',self.runname)

        self.args = mock.MagicMock()
        self.args.rundir = self.rundir
        self.args.ngsdata = self.ngsdata
        self.config = {
            'NGSDATA': self.ngsdata,
            'ion_sync': {
                'ngsdata': {
                    'default': self.ngsdata,
                    'help': 'help msg'
                }
            }
        }

    def _check_rawdata(self, rawdatapath):
        self.assertEqual(
            os.stat(self.ion_params).st_size,
            os.stat(join(rawdatapath,'ion_params_00.json')).st_size
        )
        self.assertEqual(
            sorted([basename(fq) for fq in self.fastqs]),
            sorted(os.listdir(join(rawdatapath,'plugin_out','downloads')))
        )

    def _check_readdata(self, readdatapath):
        for origpath, newname in self.samplefilemap.items():
            sympath = join(
                '../../../RawData/IonTorrent',
                self.runname,
                origpath
            )
            readpath = join(readdatapath, newname)
            self.assertTrue(
                exists(readdatapath),
                '{0} was not created correctly'.format(readdatapath)
            )
            self.assertEqual(
                sympath,
                os.readlink(readpath)
            )

    def _check_readsbysample(self, rbspath):
        for origpath, newfile in self.samplefilemap.items():
            sn, bc, rn, ext = newfile.split('.')
            sympath = join(
                '../../ReadData/IonTorrent',
                self.runname,
                newfile
            )
            readpath = join(rbspath,sn,newfile)
            self.assertTrue(
                exists(rbspath),
                '{0} does not exist or broken symlink'.format(
                    rbspath
                )
            )
            self.assertEqual(
                sympath,
                os.readlink(readpath)
            )

    def _print_ngsdata(self, ngsdatapath):
        for root, dirs, files in os.walk(ngsdatapath):
            for f in files:
                print join(root,f)

    def test_syncs_and_creates_directories(self, margparse, mconfig):
        sysargs = [self.args.rundir]
        mconfig.get_config_argparse.return_value = (
            mock.Mock(), sysargs, self.config, 'config.yaml'
        )
        margparse.ArgumentParser.return_value.parse_args.return_value = self.args
        ion_sync.main()
        self._print_ngsdata(self.ngsdata)
        self._check_rawdata(self.rawdata)
        self._check_readdata(self.readdata)
        self._check_readsbysample(self.readsbysample)

    def test_resumes_sync_and_does_not_overwrite_existing_symlink(self, margparse, mconfig):
        sysargs = [self.args.rundir]
        mconfig.get_config_argparse.return_value = (
            mock.Mock(), sysargs, self.config, 'config.yaml'
        )
        margparse.ArgumentParser.return_value.parse_args.return_value = self.args
        os.makedirs(join(self.rbs, 'sample1'))
        os.makedirs(join(self.rbs, 'sample2'))
        os.makedirs(self.readdata)
        os.makedirs(self.rawdata)
        self.rundir, self.downloads, self.ion_params = \
            self._create_rundir(self.fastqs, self.ionparam, self.rawdata)
        # This one will have readdata completed
        fqname1 = self.fastqs[0]
        src = join(
            '../../../RawData/IonTorrent',self.runname,fqname1
        )
        dst = join(
            self.readdata,self.samplefilemap[fqname1]
        )
        os.symlink(src, dst)

        # This one readdata and readsbysample are done
        fqname2 = self.fastqs[1]
        fqfile = self.samplefilemap[fqname2]
        sn, bc, rn, ext = fqfile.split('.')
        src = join(
            '../../../RawData/IonTorrent',self.runname,fqname2
        )
        src2 = join(
            '../../ReadData/IonTorrent',self.runname,fqfile
        )
        dst = join(
            self.readdata,self.samplefilemap[fqname2]
        )
        dst2 = join(
            self.rbs,sn,fqfile
        )
        os.symlink(src, dst)
        os.symlink(src2, dst2)
        with mock.patch.object(ion_sync, 'sys') as msys:
            msys.stderr = StringIO()
            ion_sync.main()
            serr = msys.stderr.getvalue().splitlines()
            print serr
            self.assertIn(
                'You may need to remove this directory and rerun the sync',
                serr[1]
            )
            self.assertIn('sample', serr[2])
            self.assertIn('ReadData', serr[2])
            self.assertIn('sample', serr[3])
            self.assertIn('ReadData', serr[3])
            self.assertIn('sample2', serr[4])
            self.assertIn('ReadsBySample', serr[4])
        self._print_ngsdata(self.ngsdata)
        self._check_rawdata(self.rawdata)
        self._check_readdata(self.readdata)
        self._check_readsbysample(self.readsbysample)

    def test_loads_defaults_from_config(self, margparse, mconfig):
        sysargs = [self.args.rundir]
        mconfig.get_config_argparse.return_value = (
            mock.Mock(), sysargs, self.config, 'config.yaml'
        )
        margparse.ArgumentParser.return_value.parse_args.return_value = self.args
        ion_sync.main()
        self._print_ngsdata(self.ngsdata)
        self._check_rawdata(self.rawdata)
        self._check_readdata(self.readdata)
        self._check_readsbysample(self.readsbysample)
