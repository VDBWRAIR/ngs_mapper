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
from nose.plugins.attrib import attr

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
            'expName': 'R_run_name',
            'experimentAnalysisSettings': {
                'barcodedSamples': json.dumps(self.barcoded_samples)
            }
        }

        # mock contents of basecaller_results
        self.basecaller_results = [
            'basecaller_results/IonXpress_001_rawlib.basecaller.bam',
            'basecaller_results/IonXpress_002_rawlib.basecaller.bam',
            'basecaller_results/IonXpress_003_rawlib.basecaller.bam',
            'basecaller_results/IonXpress_004_rawlib.basecaller.bam',
            'basecaller_results/nomatch_rawlib.basecaller.bam',
        ]

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
            'plugin_out/downloads/IonXpress_002.R_run_name.fastq':
                'sample2.IonXpress_002.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_003.R_run_name.fastq':
                'sample2.IonXpress_003.R_run_name.fastq',
            'plugin_out/downloads/IonXpress_004.R_run_name.fastq':
                'IonXpress_004.IonXpress_004.R_run_name.fastq',
            'plugin_out/downloads/nomatch.R_run_name.fastq':
                'nomatch.nomatch.R_run_name.fastq',
        }

@attr('current')
@mock.patch('__builtin__.open')
@mock.patch.object(ion_sync, 'bam')
class TestConvertBasecallerResultsToFastq(Base):
    def test_ensures_fastqoutpath_exists(self, *args):
        mapping = {
            'basecaller_results/IonXpress_001_rawlib.basecaller.bam':
                'IonXpress_001.R_run_name.fastq'
        }
        fastqpath = 'plugin_out/downloads'
        with mock.patch.object(ion_sync, 'os') as mos:
            mos.path.exists.return_value = False
            ion_sync.convert_basecaller_results_to_fastq(mapping, fastqpath)
            mos.makedirs.assert_called_once_with(fastqpath)

    def test_converts_bams_to_fastq(self, *args):
        mock_bam = args[0]
        mopen = args[1]
        fqread = '@read1\nATGC\n+\n!!!!'
        mock_bam.bam_to_fastq.return_value = iter([
            fqread
        ])
        bams = [
            'basecaller_results/IonXpress_001_rawlib.basecaller.bam',
            'basecaller_results/IonXpress_002_rawlib.basecaller.bam',
        ]
        fastqs = [
            'IonXpress_001.R_run_name.fastq',
            'IonXpress_002.R_run_name.fastq'
        ]
        mapping = dict(zip(bams,fastqs))

        fastqpath = 'plugin_out/downloads'
        ion_sync.convert_basecaller_results_to_fastq(mapping, fastqpath)
        for fq in fastqs:
            fqpath = join(fastqpath, fq)
            mopen.assert_has_call(fqpath, 'w')
            mopen.return_value.__enter__.return_value.write.assert_called_with(
                fqread + '\n'
            )

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

@attr('current')
class TestGetSamplefileMapping(Base):
    def test_creates_correct_mapping_fastq(self):
        expect = self.samplefilemap
        r = ion_sync.get_samplefile_mapping(
            self.samplemap,
            self.fastqs,
            'R_run_name'
        )
        for e in expect:
            self.assertEqual(
                expect[e],
                r[e]
            )

    def test_creates_correct_mapping_rawlibbam(self):
        expect = {
            'basecaller_results/IonXpress_001_rawlib.basecaller.bam':
                'sample1.IonXpress_001.R_run_name.fastq',
            'basecaller_results/IonXpress_004_rawlib.basecaller.bam':
                'IonXpress_004.IonXpress_004.R_run_name.fastq'
        }
        r = ion_sync.get_samplefile_mapping(
            self.samplemap,
            expect.keys(),
            'R_run_name'
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
            self.fastqs,
            'R_run_name'
        )

@attr('current')
@mock.patch.object(ion_sync, 'glob')
@mock.patch.object(ion_sync, 'json')
@mock.patch('__builtin__.open', mock.Mock())
class TestIonMapping(Base):
    def test_ion_params_maps_pluginout_downloads_fastq(self, mjson, mglob):
        mjson.load.return_value = self.ionparam
        mjson.loads.return_value = self.barcoded_samples
        mglob.glob.return_value = self.fastqs
        r = ion_sync.ion_mapping(
            self.fastqs,
            'ion_param_00.json'
        )
        for origpath in self.samplefilemap:
            self.assertEqual(
                self.samplefilemap[origpath],
                r[origpath]
            )

    def test_ion_params_maps_basecaller_bams(self, mjson, mglob):
        mjson.load.return_value = self.ionparam
        mjson.loads.return_value = self.barcoded_samples
        mglob.glob.return_value = self.basecaller_results
        samplefilemap = {
            'basecaller_results/IonXpress_001_rawlib.basecaller.bam':
                'sample1.IonXpress_001.R_run_name.fastq',
            'basecaller_results/IonXpress_002_rawlib.basecaller.bam':
                'sample2.IonXpress_002.R_run_name.fastq'
        }
        r = ion_sync.ion_mapping(
            samplefilemap.keys(),
            'ion_param__00.json'
        )
        for origpath in samplefilemap:
            self.assertEqual(
                samplefilemap[origpath],
                r[origpath]
            )

    def test_ion_params_missing_expName(self, mjson, mglob):
        del self.ionparam['expName']
        self.assertRaises(
            ion_sync.InvalidIonParam,
            self.test_ion_params_maps_pluginout_downloads_fastq,
        )
        self.ionparam['expName'] = ''
        self.assertRaises(
            ion_sync.InvalidIonParam,
            self.test_ion_params_maps_pluginout_downloads_fastq,
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
        if path is none than make tempdir

        makes:
            path/basecaller_results/ with reads in it
            path/plugin_out/downloads with reads in it if reads has that path
            path/ion_param__00.json

        reads is a list of 
        plugin_out/downloads/IonXpress_XXX.whatever.fastq
            or
        basecaller_results/IonXpress_XXX_rawlib.basecaller.bam

        ionparam is a dict representing the json.load from actual
        ion_params_00.json file
        '''
        if path is None:
            rundir = tempfile.mkdtemp(prefix='rundir_',suffix='ion',dir=self.tdir)
        else:
            rundir = path
        downloads = join(rundir, dirname(reads[0]))
        ion_params = join(rundir,'ion_params_00.json')
        basecaller_results_dir = join(rundir, 'basecaller_results')
        # Make plugin_out/downloads dir
        os.makedirs(downloads)
        # If reads is just the bam files location then it may be already created
        if not isdir(basecaller_results_dir):
            os.makedirs(basecaller_results_dir)
            self._make_basecaller_results(rundir, reads)
        # Make ion_params_00.json file
        with open(ion_params,'w') as fh:
            json.dump(ionparam, fh)
        self._make_plugin_out_downloads(rundir, reads)
        return rundir, downloads, ion_params

    def _make_basecaller_results(self, outdir, reads):
        ''' make basecaller_results regardless of dirname of reads '''
        for read in reads:
            # Just ionxpress part
            m = re.search('(IonXpress_\d+|rawtf|nomatch)', read).group(1)
            bam = join(outdir, 'basecaller_results', m + '_rawlib.basecaller.bam')
            open(bam,'w').close()

    def _make_plugin_out_downloads(self, outdir, reads):
        # Make mock fastqs
        for read in reads:
            readpath = join(outdir,read)
            open(readpath,'w').close()

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
            '/missing/path/IonXpress_001.broken.bam',
            join(self.downloads,'IonXpress_001.bam')
        )
        self.runname = basename(self.rundir)
        self.readdata = join(self.readdata,'IonTorrent',self.runname)
        self.rawdata = join(self.rawdata,'IonTorrent',self.runname)

        self.args = mock.MagicMock()
        self.args.rundir = self.rundir
        self.args.ngsdata = self.ngsdata
        self.args.print_samplemapping = False
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
        e = sorted([basename(fq) for fq in self.fastqs])
        print e
        r = sorted(os.listdir(join(rawdatapath,'plugin_out','downloads')))
        print r
        self.assertEqual(e, r)

    def _check_readdata(self, readdatapath, samplefilemap):
        for origpath, newname in samplefilemap.items():
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

    def _print_path(self, ngsdatapath):
        for root, dirs, files in os.walk(ngsdatapath):
            for f in files:
                print join(root,f)

    def test_only_prints_samplemapping(self, margparse, mconfig):
        sysargs = [self.args.rundir, '--print-samplemapping']
        mconfig.get_config_argparse.return_value = (
            mock.Mock(), sysargs, self.config, 'config.yaml'
        )
        self.args.print_samplemapping = True
        margparse.ArgumentParser.return_value.parse_args.return_value = self.args
        with mock.patch.object(ion_sync, 'sys') as msys:
            with mock.patch.object(ion_sync, 'sync_readdata') as msync_readdata:
                with mock.patch.object(ion_sync, 'sync_readsbysample') as msync_readsbysample:
                    with mock.patch.object(ion_sync, 'shutil') as mshutil:
                        ion_sync.main()
                        self.assertFalse(msync_readdata.called)
                        self.assertFalse(msync_readsbysample.called)
                        self.assertFalse(mshutil.copytree.called)
                        for k,v in self.samplefilemap.items():
                            msys.stdout.write.assert_has_call("{0} -> {1}\n".format(k,v))

    @attr('current')
    def test_converts_basecaller_results_bams_to_fastq(self, margparse, mconfig):
        # Remove fastq to make sure they are rebuilt in ReadData
        plugin_out_dir = join(self.rundir, 'plugin_out', 'downloads')
        for f in glob(join(plugin_out_dir,'*')):
            os.unlink(f)

        sysargs = [self.args.rundir]
        mconfig.get_config_argparse.return_value = (
            mock.Mock(), sysargs, self.config, 'config.yaml'
        )
        margparse.ArgumentParser.return_value.parse_args.return_value = self.args
        ion_sync.main()

        self._print_path(self.ngsdata)
        print '-----'
        self._print_path(self.tdir)

        self._check_rawdata(self.rawdata)
        self._check_readdata(self.readdata, self.samplefilemap)
        self._check_readsbysample(self.readsbysample)

    def test_syncs_and_creates_directories(self, margparse, mconfig):
        sysargs = [self.args.rundir]
        mconfig.get_config_argparse.return_value = (
            mock.Mock(), sysargs, self.config, 'config.yaml'
        )
        margparse.ArgumentParser.return_value.parse_args.return_value = self.args
        ion_sync.main()
        self._check_rawdata(self.rawdata)
        self._check_readdata(self.readdata, self.samplefilemap)
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
        self._print_path(self.ngsdata)
        self._check_rawdata(self.rawdata)
        self._check_readdata(self.readdata, self.samplefilemap)
        self._check_readsbysample(self.readsbysample)

    def test_loads_defaults_from_config(self, margparse, mconfig):
        sysargs = [self.args.rundir]
        mconfig.get_config_argparse.return_value = (
            mock.Mock(), sysargs, self.config, 'config.yaml'
        )
        margparse.ArgumentParser.return_value.parse_args.return_value = self.args
        ion_sync.main()
        self._print_path(self.ngsdata)
        self._check_rawdata(self.rawdata)
        self._check_readdata(self.readdata, self.samplefilemap)
        self._check_readsbysample(self.readsbysample)
