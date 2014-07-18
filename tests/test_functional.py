from imports import *
from test_install import Base
import json

# Derive from installer base to make life easier
class BaseFunctional(Base):
    @classmethod
    def setUpClass( klass ):
        super(Base,klass).setUpClass()
        # Uninstall first since __init__.py already installed for us
        # Only install once because it takes a long time
        #klass.returncode, klass.output = klass.run_installer( klass.vpath )

    @classmethod
    def compile_functional_fixtures( klass, fixdir ):
        ''' Return list of [(readdir,conf),...] '''
        fixtures = []
        for df in os.listdir( fixdir ):
            df = join( fixdir, df )
            if isdir(df):
                confpath = join( dirname(df), basename(df) + '.conf' )
                fixtures.append( (df, confpath) )
        return fixtures

    @classmethod
    def parse_conf( klass, confpath ):
        conf = {}
        with open(confpath) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith('#'):
                    continue
                k,v = line.split(':', 1)
                conf[k] = v
        return conf

    @classmethod
    def make_samplesheet( klass, fixtures, outputpath ):
        with open(outputpath,'w') as fh:
            for readsdir, conf in fixtures:
                c = klass.parse_conf( conf )
                ref = join( dirname(readsdir), c['reference'] )
                fh.write( '{}\t{}\n'.format(basename(readsdir),ref) )

    @classmethod
    def run_fixtures( klass, fixtures ):
        samplesheet = 'ss.tsv'
        klass.make_samplesheet(fixtures,'ss.tsv')
        runsh_sh = TestRunPipeline.script_path( 'runsamplesheet.sh' )
        # All fixtures should be in same dir, so just grab the dirname of the first
        rbsdir = dirname( fixtures[0][0] )
        cmd = '{} {} {}'.format(runsh_sh,rbsdir,samplesheet)
        ret,out = TestRunPipeline.run_script( cmd )

        return samplesheet,ret,out

class TestRunPipeline(BaseFunctional):
    @classmethod
    def setUpClass( klass ):
        klass.fixtures = klass.compile_functional_fixtures( join(fixtures.FIXDIR,'functional') )
        klass.samplesheet, klass.ret, klass.out = klass.run_fixtures( klass.fixtures )
        #klass.samplesheet, klass.ret, klass.out = 'ss.tsv',0,''

    def setUp( self ):
        self.fixtures = self.__class__.fixtures
        self.samplesheet = self.__class__.samplesheet
        self.returncode = self.__class__.ret
        self.output = self.__class__.out

    def test_return_code_and_output(self):
        eq_( 0, self.returncode, 'Return code from running runsamplesheet.sh was not 0' )
        for reads, config in self.fixtures:
            sn = basename(reads)
            p = 'Please check the logfile (/dev/shm/\w+/{}.log)'.format(sn)
            m = re.search( p, self.output, re.S|re.M )
            if m:
                print m.group(1)
                print open(m.group(1)).read()
            ok_( 'Starting {}'.format(sn) in self.output, "Did not start {}".format(sn) )
            ok_( 'Finished {}'.format(sn) in self.output, "Did not finish {}".format(sn) )

    def check_sample_project_files( self, projdir, fixture ):
        # Files defined that should exist
        efiles = self.__class__.parse_conf( fixture[1] )['files'].split()
        failed = [join(projdir,ef) for ef in efiles if not exists(join(projdir,ef))]
        for f in failed:
            print 'Pipeline did not produce project file {}'.format(f)
        if failed:
            ok_( False )

    def test_pipeline_produced_expected_files_dirs( self ):
        f = 'file'
        d = 'directory'
        expected = [
            ('graphsample.log',f),
            ('MapUnmapReads.png',f),
            ('pipeline.log',f),
            ('PipelineTimes.png',f),
            ('QualDepth.pdf',f),
            ('ss.tsv',f),
            ('Projects',d),
            ('vcf_consensus',d),
        ]

        for e, typ in expected:
            ok_( exists( e ), 'Pipeline did not produce {}'.format(e) )

    def test_project_directories_have_expected_files( self ):
        # Ensure each project has correct files too
        for fixture in self.fixtures:
            sn = basename(fixture[0])
            projdir = join( 'Projects', sn )
            self.check_sample_project_files( projdir, fixture )

    def test_vcf_consensus_has_symlink_consensuses( self ):
        for pdir in glob( join('Projects','*') ):
            sn = basename(pdir)
            consensus_file = join( pdir, sn + '.bam.consensus.fasta' )
            consensus_link = join( 'vcf_consensus', sn + '.fasta' )
            ok_( samefile(consensus_file, consensus_link), 'vcf_consensus file {} is not correctly linked'.format(consensus_link) )

    def count_mutations( self, fastapath ):
        amb_bases = 'nmrwsykvhdb'
        counts = {}
        for seq in SeqIO.parse(fastapath,'fasta'):
            id = seq.description
            # Init counts to 0 for this reference
            counts[id] = {b:0 for b in amb_bases}
            # Count amb bases
            for b in seq.seq:
                b = b.lower()
                if b in amb_bases:
                    counts[id][b] += 1
        return counts

    def get_fixture_mutation_counts( self, fixture ):
        c = self.__class__.parse_conf( fixture[1] )
        counts = c['mutations']
        counts = json.loads( counts )
        return counts

    def test_consensus_mutations( self ):
        for reads,c in self.fixtures:
            sn = basename(reads)
            consensus_file = join( 'vcf_consensus', sn + '.fasta' )
            # Count all mutations in resulting project
            mutation_counts = self.count_mutations( consensus_file )
            # Retrieve the expected mutation count
            expected_counts = self.get_fixture_mutation_counts( (reads,c) )

            # Compare the two counts
            for ref, counts in expected_counts.iteritems():
                result_counts = mutation_counts[ref]
                for amb_base, count in counts.iteritems():
                    eq_( count, result_counts[amb_base] )
