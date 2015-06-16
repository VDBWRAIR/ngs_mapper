import tempfile
import shutil
import os
from os.path import *
from . import tdir
import subprocess

class BaseTester(object):
    def _C( self, *args, **kwargs ):
        '''
        Set modulepath as instance variable to be the name of the module
        Set functionname as instance variable to automagically run that function
        with self._C
        '''
        m = __import__( self.modulepath, fromlist=[self.functionname] )
        return getattr(m,self.functionname)( *args, **kwargs )

class BaseClass( BaseTester ):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        self.tempdir = tempfile.mkdtemp(prefix='unit',suffix='test',dir=tdir)
        os.chdir(self.tempdir)
        self.mpileups = {
            'Ref1': [
                self._mock_pileup_str('Ref1', 1, 'A', 1, 'A', '!!', '!!'),
            ],
            'Ref2': [
                self._mock_pileup_str('Ref2', 1, 'A', 2, 'CC', '!!', '!!'),
                self._mock_pileup_str('Ref2', 2, 'A', 2, 'CC', '!!', '!!'),
            ],
            'Ref3': [
                self._mock_pileup_str('Ref3', 1, 'A', 3, 'TTT', '!!!', '!!!'),
                self._mock_pileup_str('Ref3', 2, 'A', 3, 'TTT', '!!!', '!!!'),
                self._mock_pileup_str('Ref3', 3, 'A', 3, 'TTT', '!!!', '!!!'),
            ],
            'Ref4': [
            ],
            'Ref5': [
                self._mock_pileup_str('Ref5', 1, 'A', 1, 'G', '!', '!'),
                self._mock_pileup_str('Ref5', 3, 'A', 1, 'G', '!', '!'),
                self._mock_pileup_str('Ref5', 4, 'A', 1, 'G', '!', '!'),
                self._mock_pileup_str('Ref5', 5, 'A', 1, 'G', '!', '!'),
            ]
        }

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
        return script

    @classmethod
    def run_script( self, script ):
        print "Running {0}".format(script)
        try:
            return (0,subprocess.check_output( script, stderr=subprocess.STDOUT, shell=True ))
        except subprocess.CalledProcessError as e:
            return (e.returncode, e.output)

    def _mock_pileup_str(self, *args):
        return '\t'.join([str(x) for x in args])

import fixtures
from ngs_mapper.bam import indexbam
class BaseBamRef(BaseClass):
    bam = join(fixtures.THIS,'fixtures','varcaller','paired.bam.gz')
    ref = join(fixtures.THIS,'fixtures','varcaller','ref.fasta.gz')
    mytempdir = ''

    @classmethod
    def setUpClass(klass):
        # Unpacks everything once so it doesn't slow down so much
        super(BaseBamRef,klass).setUpClass()
        import tempfile
        klass.mytempdir = tempfile.mkdtemp(prefix='basebamref',suffix='test',dir=tdir)
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

def make_seqrec( seq, quals, id='id' ):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    seq = Seq( seq, generic_dna )
    rec = SeqRecord(
        seq,
        id=id,
        description='description',
        name='name'
    )
    rec._per_letter_annotations['phred_quality'] = quals
    return rec

def random_seqs( numseqs=100 ):
    import random
    dna = 'ATGC'
    seqs = []
    maxlen = 0.0
    maxqual = 0.0
    for i in range( 1, numseqs ):
        randbase = random.choice( dna )
        randlen = random.randint( 0, 1000 )
        randseq = ''.join( [random.choice(dna) for i in range(0,randlen)] )
        randqual = [random.randint(0,60) for i in range(0, randlen)]
        aqual = 0
        if randlen != 0:
            aqual = round( sum(randqual) * 1.0 / randlen )
            maxqual = max( aqual, maxqual )
        maxlen = max( maxlen, randlen )
        seqs.append( make_seqrec( randseq, randqual ) )
    return (seqs, maxlen, maxqual)

def rand_seqrec( seqlen, cal, car, cql, cqr ):
    '''
        Makes a Bio.SeqRecord.SeqRecord such that the .seq returns
        a Bio.Seq.Seq with a random sequence that has length seqlen if seqlen is int otherwise makes Seq with seqlen
        Then the SeqRecord has ._per_letter_annotations['phred_quality'] with random qualities
        for all the bases
        SeqRecord.annotations clip_adapter_left, right, clip_qual_left, right are set to cal, car, cql and cqr
    '''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_dna
    import random
    seq = None
    if isinstance( seqlen, int ):
        # Random Sequence
        seq = Seq( rand_seq(seqlen), generic_dna )
        # Random qualities
        qual = [random.randint(1,40) for i in range(seqlen)]
    else:
        seq = Seq( seqlen, generic_dna )
        # Random qualities
        qual = [random.randint(1,40) for i in range(len(seqlen))]
    # Random id. Hopefully random enough so no duplicate ids
    id = 'seq_{0}'.format(random.randint(1,999999999999))
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

def touch( filepath ):
    subprocess.call(['touch',filepath])
