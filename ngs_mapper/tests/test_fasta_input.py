import sh
import unittest
import tempfile
from os.path import join, dirname, abspath
import shutil
from Bio import SeqIO
from functools import partial
import plumbum
THISD = dirname(abspath(__file__))
here = partial(join, THISD)
def runsample(indir, outdir):

    plumbum.local['runsample'][indir, here("fixtures/functional/947.ref.fasta"), \
                 "947", "-od", outdir, "--fasta"]()
#    sh.runsample(indir, here("fixtures/functional/947.ref.fasta"), \
#                 "947", index_min=0, outdir=outdir, fasta=True)

class TestFastaInput(unittest.TestCase):

    def setUp(self):
        self.fastaInputDir = tempfile.mkdtemp()
        self.fastqInputDir = tempfile.mkdtemp()
        self.fastaOutputDir = join(dirname(self.fastaInputDir), 'fastaout')
        self.fastqOutputDir = join(dirname(self.fastqInputDir), 'fastqout')
        fq = here("fixtures/functional/947/947_S32_L001_R1_001_2013_12_17.fastq")
        fa = "R1.fasta"
        shutil.copy(fq, self.fastqInputDir)
        SeqIO.convert(fq, 'fastq', join(self.fastaInputDir, fa), 'fasta')

    def tearDown(self):
       shutil.rmtree(self.fastaInputDir )
       shutil.rmtree(self.fastqInputDir )
       shutil.rmtree(self.fastaOutputDir)
       shutil.rmtree(self.fastqOutputDir)
       pass

    def test_run_stats_with_fasta_equals_run_with_fastq_40s(self):
        runsample(self.fastaInputDir, self.fastaOutputDir)
        runsample(self.fastqInputDir, self.fastqOutputDir)
        fastaStats = join(self.fastaOutputDir, "flagstats.txt")
        fastqStats = join(self.fastqOutputDir, "flagstats.txt")
        fastaStats = open(fastaStats).read()
        fastqStats = open(fastqStats).read()
        self.assertEquals(fastqStats, fastaStats)
