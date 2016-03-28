import sh
import unittest
import tempfile
from os.path import join, dirname, abspath
import shutil
from Bio import SeqIO
from functools import partial
THISD = dirname(abspath(__file__))
here = partial(join, THISD)
inputFastq = "fixtures/fasta/780_S12_L001_R1_001_2014_04_16.fastq"
def runsample(indir, outdir):
    sh.runsample(indir, here("fixtures/functional/780.ref.fasta"), \
                 "780", outdir=outdir, fasta=True)

class TestFastaInput(unittest.TestCase):

    def setUp(self):
        self.fastaInputDir = tempfile.mkdtemp()
        self.fastqInputDir = tempfile.mkdtemp()
        self.fastaOutputDir = join(dirname(self.fastaInputDir), 'fastaout')
        self.fastqOutputDir = join(dirname(self.fastqInputDir), 'fastqout')
        fq = here(inputFastq)
        fa = "R1.fasta"
        shutil.copy(fq, self.fastqInputDir)
        SeqIO.convert(fq, 'fastq', join(self.fastaInputDir, fa), 'fasta')

    def tearDown(self):
       shutil.rmtree(self.fastaInputDir )
       shutil.rmtree(self.fastqInputDir )
#       shutil.rmtree(self.fastaOutputDir)
#       shutil.rmtree(self.fastqOutputDir)
       pass

    def test_run_stats_with_fasta_equals_run_with_fastq_40s(self):
        runsample(self.fastaInputDir, self.fastaOutputDir)
        runsample(self.fastqInputDir, self.fastqOutputDir)
        fastaStats = join(self.fastaOutputDir, "flagstats.txt")
        fastqStats = join(self.fastqOutputDir, "flagstats.txt")
        fastaStats = open(fastaStats).read()
        fastqStats = open(fastqStats).read()
        self.assertEquals(fastqStats, fastaStats)
