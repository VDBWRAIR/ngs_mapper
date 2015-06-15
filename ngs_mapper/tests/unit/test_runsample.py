from . import unittest
import mock

from ngs_mapper import runsample

class TestPBSJob(unittest.TestCase):
    def setUp(self):
        self.runsample_args = self.rsa = [
            'readsdir', 'foo.fasta', 'prefix',
            '-trim_qual', '25', '-head_crop', '0',
            '-minth', '0', '--CN', 'bar', '-od', 'outdir'
        ]
        self.qsub_args = self.qa = [
            '--qsub_M', 'email@example.com',
            '--qsub_l', 'nodes=1:ppn=1'
        ]
        self.args = [
            self.runsample_args,
            self.qsub_args
        ]

    def test_calls_runsample_correctly(self):
        r = runsample.pbs_job(*self.args)
        self.assertIn(
            'cd $PBS_O_WORKDIR',
            r
        )
        self.assertIn(
            'runsample readsdir foo.fasta prefix -trim_qual 25 -head_crop 0 ' \
            '-minth 0 --CN bar -od outdir',
            r
        )
        self.assertNotIn(
            'export TMPDIR',
            r
        )

    def test_sets_correct_pbs_directives(self):
        r = runsample.pbs_job(*self.args)
        self.assertIn(
            '#PBS -N prefix-ngs_mapper',
            r
        )
        self.assertIn(
            '#PBS -M email@example.com',
            r
        )
        self.assertIn(
            '#PBS -l nodes=1:ppn=1',
            r
        )

    def test_omits_email_if_not_in_args(self):
        self.args[1] = ['--qsub_l', 'nodes=1:ppn=1']
        r = runsample.pbs_job(*self.args)
        self.assertNotIn('#PBS -M', r)
        self.assertNotIn('#PBS -m', r)

    def test_sets_TMPDIR_only_if_already_in_environment(self):
        with mock.patch.dict(runsample.os.environ, TMPDIR='/path/foo'):
            r = runsample.pbs_job(*self.args)
            self.assertIn(
                'export TMPDIR=/path/foo',
                r
            )

class TestSplitArgs(unittest.TestCase):
    def test_splits_correctly(self):
        r = runsample.split_args(
            'readsdir refpath prefix --qsub_l bla=1:ble=5 -foo bar -baz jaz ' \
            '--qsub_M user@example.com'
        )
        self.assertEqual(
            'readsdir refpath prefix -foo bar -baz jaz'.split(),
            r[0]
        )
        self.assertEqual(
            '--qsub_l bla=1:ble=5 --qsub_M user@example.com'.split(),
            r[1]
        )

    def test_returns_empty_qsub_args_if_none(self):
        r = runsample.split_args(
            'readsdir refpath prefix -foo bar -baz jaz'
        )
        self.assertEqual(
            'readsdir refpath prefix -foo bar -baz jaz'.split(),
            r[0]
        )
        self.assertEqual([], r[1])
