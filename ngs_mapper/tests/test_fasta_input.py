import sh
import unittest
from glob import glob
import tempfile
from os.path import join, dirname
import shutil
from Bio import SeqIO

def runsample(indir, outdir):
    sh.runsample(indir, "tests/fixtures/functional/947.ref.fasta", "947", od=outdir)

class TestFastaInput(unittest.TestCase):

    def setUp(self):
        self.fastaInputDir = tempfile.mkdtemp()
        self.fastqInputDir = tempfile.mkdtemp()
        self.fastaOutputDir = join(dirname(self.fastaInputDir), 'fastaout')
        self.fastqOutputDir = join(dirname(self.fastqInputDir), 'fastqout')
        fq = "tests/fixtures/functional/947/947_S32_L001_R1_001_2013_12_17.fastq"
        fa = "R1.fasta"
        shutil.copy(fq, self.fastqInputDir)
        SeqIO.convert(fq, 'fastq', join(self.fastaInputDir, fa), 'fasta')

    def tearDown(self):
       shutil.rmtree(self.fastaInputDir )
       shutil.rmtree(self.fastqInputDir )
       shutil.rmtree(self.fastaOutputDir)
       shutil.rmtree(self.fastqOutputDir)

    def test_run_stats_with_fasta_equals_run_with_fastq_40s(self):
        runsample(self.fastaInputDir, self.fastaOutputDir)
        runsample(self.fastqInputDir, self.fastaOutputDir)
        fastaStats = join(self.fastaOutputDir, "flagstats.txt")
        fastqStats = join(self.fastaOutputDir, "flagstats.txt")
        fastaStats = open(fastaStats).read()
        fastqStats = open(fastqStats).read()
        self.assertEquals(fastqStats, fastaStats)
