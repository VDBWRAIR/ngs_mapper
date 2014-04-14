import tempfile
import shutil
import os
from os.path import *
from . import tdir
import subprocess

class BaseClass( object ):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        self.tempdir = tempfile.mkdtemp(prefix='unit',suffix='test',dir=tdir)
        os.chdir(self.tempdir)

    def tearDown(self):
        os.chdir(tdir)

    def create_file(self,filepath,contents):
        linecount = 0
        with open(filepath,'w') as fh:
            for line in contents.splitlines(True):
                fh.write(line)
                linecount += 1
        return linecount

    @classmethod
    def script_path( self, script ):
        path = dirname( dirname( __file__ ) )
        return join( path, script )

    @classmethod
    def run_script( self, script ):
        print "Running {}".format(script)
        try:
            return (0,subprocess.check_output( script, stderr=subprocess.STDOUT, shell=True ))
        except subprocess.CalledProcessError as e:
            return (e.returncode, e.output)

import fixtures
from bam import indexbam
class BaseBamRef(BaseClass):
    bam = join(fixtures.THIS,'fixtures','varcaller','paired.bam.gz')
    ref = join(fixtures.THIS,'fixtures','varcaller','ref.fasta.gz')
    mytempdir = ''

    @classmethod
    def setUpClass(klass):
        # Unpacks everything once so it doesn't slow down so much
        super(BaseBamRef,klass).setUpClass()
        import tempfile
        klass.mytempdir = tempfile.mkdtemp(prefix='basebamref',suffix='test')
        klass.bam = fixtures.ungiz(klass.bam,klass.mytempdir)
        klass.ref = fixtures.ungiz(klass.ref,klass.mytempdir)
        klass.bamindex = indexbam( klass.bam )

    @classmethod
    def tearDownClass(klass):
        super(BaseBamRef,klass).tearDownClass()
        import shutil
        shutil.rmtree(klass.mytempdir)

class BaseBaseCaller(BaseClass):
    def setUp( self ):
        super( BaseBaseCaller, self ).setUp()
        fixpath = join( fixtures.THIS, 'fixtures', 'base_caller' )
        self.bam = join( fixpath, 'test.bam' )
        self.bai = join( fixpath, 'test.bam.bai' )
        self.ref = join( fixpath, 'testref.fasta' )
        self.sam = join( fixpath, 'test.sam' )
        self.vcf = join( fixpath, 'test.vcf' )
        self.template = join( fixpath, 'template.vcf' )

def rand_seqrec( seqlen, cal, car, cql, cqr ):
    '''
        Makes a random Bio.SeqRecord.SeqRecord such that the .seq returns
        a Bio.Seq.Seq with a random sequence that has length seqlen
        Then the SeqRecord has ._per_letter_annotations['phred_quality'] with random qualities
        for all the bases
        SeqRecord.annotations clip_adapter_left, right, clip_qual_left, right are set to cal, car, cql and cqr
    '''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    import random
    # Random Sequence
    seq = Seq( rand_seq(seqlen), generic_dna )
    # Random qualities
    qual = [random.randint(1,40) for i in range(seqlen)]
    # Random id. Hopefully random enough so no duplicate ids
    id = 'seq_{}'.format(random.randint(1,999999999999))
    record = SeqRecord(
        seq,
        id=id,
        description='Random sequence',
        name=id
    )
    record._per_letter_annotations['phred_quality'] = qual
    record.annotations['clip_adapter_left'] = cal
    record.annotations['clip_adapter_right'] = car
    record.annotations['clip_qual_left'] = cql
    record.annotations['clip_qual_right'] = cqr
    return record

def rand_seq( seqlen ):
    ''' return random sequence length seqlen '''
    import random
    dna = ('A','C','G','T')
    return ''.join( [dna[random.randint(0,3)] for i in range(seqlen)] )
