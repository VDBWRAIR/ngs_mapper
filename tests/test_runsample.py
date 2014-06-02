from imports import *

class Base(common.BaseBamRef):
    reads_by_sample = ''
    f = ''
    r = ''
    @classmethod
    def setUpClass(klass):
        # Sets up a ReadsBySample directory and a reference to map to
        super(Base,klass).setUpClass()
        # Setup readsbysample directory
        klass.reads_by_sample = join(THIS,'fixtures','reads')

        # Fetch the reference
        _,__,ref = fixtures.get_sample_paired_reads_compressed()

        # Rename ref to correct extension
        ref = fixtures.ungiz( ref, klass.mytempdir )
        klass.ref = join( dirname(ref), 'reference.fa' )
        os.rename( ref, klass.ref )

    def check_git_repo( self, path ):
        gitdir = join( path, '.git' )
        cmd = 'git status status'.format( path, gitdir )
        try:
            output = subprocess.check_output( cmd, cwd=path, stderr=subprocess.STDOUT, shell=True )
            return 'Not a git repository' not in output
        except subprocess.CalledProcessError as e:
            print e.output
            return False

class TestUnitArgs(object):
    def _C( self, arglist ):
        from runsample import parse_args
        return parse_args( arglist )

    def test_defaults( self ):
        args = ['ReadsBySample','Reference.fasta','Sample1']
        res = self._C( args )
        eq_( 'ReadsBySample', res.readsdir )
        eq_( 'Reference.fasta', res.reference )
        eq_( 'Sample1', res.prefix )
        eq_( os.getcwd(), res.outdir )

    def test_set_outdir( self ):
        args = ['ReadsBySample','Reference.fasta','Sample1','-od','outdir']
        res = self._C( args )
        eq_( 'ReadsBySample', res.readsdir )
        eq_( 'Reference.fasta', res.reference )
        eq_( 'Sample1', res.prefix )
        eq_( 'outdir', res.outdir )

@patch('runsample.logger',Mock())
class TestUnitRunCMD(object):
    import runsample
    def _C( self, cmdstr, stdin=None, stdout=None, stderr=None, script_dir=None ):
        from runsample import run_cmd
        if script_dir is None:
            return run_cmd( cmdstr, stdin, stdout, stderr )
        else:
            return run_cmd( cmdstr, stdin, stdout, stderr, script_dir )

    def test_runs_command( self ):
        import subprocess
        res = self._C( '/bin/echo test', stdout=subprocess.PIPE )
        eq_( ('test\n',None), res.communicate() )

    def test_does_pipe( self ):
        import subprocess
        p1 = self._C( '/bin/echo test', stdout=subprocess.PIPE )
        p2 = self._C( '/bin/cat -', stdin=p1.stdout, stdout=subprocess.PIPE )
        p1.stdout.close()
        eq_( ('test\n',None), p2.communicate() )

    def test_script_dir_none( self ):
        self._C( 'echo', script_dir='' )

    def test_script_dir_somepath( self ):
        self._C( 'echo', script_dir='/bin' )

    @raises(runsample.MissingCommand)
    def test_missing_cmd_exception_caught( self ):
        self._C( 'missing.sh' )

@patch('os.path.isdir')
@patch('os.path.exists')
class TestUnitTempProjdir(object):
    def _C( self, suffix, prefix ):
        from runsample import temp_projdir
        return temp_projdir( suffix, prefix )

    def test_has_tempfs( self, exists, isdir ):
        exists.return_value = True
        isdir.return_value = True
        res = self._C( 'shmtest', 'test' )
        d, bn = split( res )
        eq_( '/dev/shm', d )
        assert bn.startswith( 'shmtest' )
        assert bn.endswith( 'test' )

    def test_nothas_tempfs( self, exists, isdir ):
        exists.return_value = False
        isdir.return_value = False
        res = self._C( 'shmtest', 'test' )
        d, bn = split( res )
        eq_( '/tmp', d )
        assert bn.startswith( 'shmtest' )
        assert bn.endswith( 'test' )

@attr('current')
@patch('runsample.logger',Mock())
class TestMakeProjectRepo(Base):
    def _C( self, *args, **kwargs ):
        from runsample import make_project_repo
        return make_project_repo( *args, **kwargs )

    def test_no_existing_repo( self ):
        path = 'outdir'
        os.mkdir( path ) 
        self._C( path )
        ok_( self.check_git_repo( path ) )
    
    def test_existing_repo( self ):
        path = 'outdir'
        os.mkdir( path ) 
        curd = os.getcwd()
        os.chdir( path )
        subprocess.call( 'git init', shell=True )
        os.chdir( curd )
        self._C( path )
        ok_( self.check_git_repo( path ) )

class TestFunctional(Base):
    def _run_runsample( self, readdir, reference, fileprefix, od=None ):
        script_path = dirname( dirname( abspath( __file__ ) ) )
        script_path = join( script_path, 'runsample.py' )
        cmd = script_path + ' {} {} {}'.format(readdir, reference, fileprefix)
        if od is not None:
            cmd += ' -od {}'.format(od)
        print "Running: {}".format(cmd)
        cmd = shlex.split( cmd )
        try:
            sout = subprocess.check_output( cmd, stderr=subprocess.STDOUT )
        except subprocess.CalledProcessError as e:
            return (e.output,-1)
        return sout,0

    def _ensure_expected_output_files( self, outdir, prefix ):
        efiles = self._expected_files( outdir, prefix )
        ef = set( [x for y,x in efiles] )
        rf = set( [join(outdir,f) for f in os.listdir( outdir )] )
        print "Files missing from project:"
        print ef - rf
        print "Extra files in project:"
        print rf - ef
        eq_( ef, rf )
        for typ, ef in efiles:     
            if typ == 'file':
                assert isfile( ef ), "{} was not created".format(ef)
                assert os.stat( ef ).st_size > 0, "{} was not > 0 bytes".format(ef)
            else:
                assert isdir( ef ), "{} was not created".format(ef)

    def _expected_files( self, outdir, prefix ):
        efiles = []
        bamfile = join( outdir, prefix + '.bam' )
        f = 'file'
        d = 'directory'
        efiles.append( (f,bamfile) )
        efiles.append( (f,bamfile + '.bai') )
        efiles.append( (f,bamfile + '.qualdepth.json') )
        efiles.append( (f,bamfile + '.qualdepth.png') )
        efiles.append( (f,bamfile + '.consensus.fasta') )
        efiles.append( (f,join( outdir, 'bwa.log') ) )
        efiles.append( (f,join( outdir, 'flagstats.txt') ) )
        efiles.append( (f,join( outdir, prefix + '.std.log') ) )
        efiles.append( (f,join( outdir, prefix + '.log') ) )
        efiles.append( (f,bamfile + '.vcf') )
        efiles.append( (d,join( outdir, 'qualdepth') ) )
        efiles.append( (d,join( outdir, 'trimmed_reads' )) )
        efiles.append( (f,join( outdir, prefix+'.reads.png' )) )
        efiles.append( (d,(join( outdir, '.git' ))) )

        # Reference and indexes
        ref = join(outdir, basename(self.ref))
        indexes = [(f,ref +'.'+ e) for e in ('sa','amb','ann','bwt','pac')]
        efiles.append( (f,(join( outdir, basename(self.ref)))) )
        efiles += indexes

        return efiles

    @attr('current')
    def test_runs_correctly( self ):
        projdir = 'outdir'
        prefix = 'testsample'
    
        out,ret = self._run_runsample( self.reads_by_sample, self.ref, prefix, projdir )
        ok_( self.check_git_repo( projdir ), 'Did not create Git repository for project' )
        self._ensure_expected_output_files( projdir, prefix )

    def test_ensure_samplename_in_consensus( self ):
        out,ret = self._run_runsample( self.reads_by_sample, self.ref, 'testsample', 'outdir' )
        from Bio import SeqIO
        cons_file = join( 'outdir', 'testsample.bam.consensus.fasta' )
        # Every seq should have same id as the prefix passed to the command
        for seq in SeqIO.parse( cons_file, 'fasta' ):
            eq_( 'testsample', seq.id )

    def test_outdir_exists_nonempty_should_skip( self ):
        os.mkdir( 'outdir' )
        res,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        eq_( 0, ret )
        assert 'AlreadyExists' not in res, "Raises exception when it should not have"
        res,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        eq_( -1, ret )
        assert 'AlreadyExists' in res, "Did not raise exception"

    def test_outdir_exists_empty( self ):
        os.mkdir( 'outdir' )
        out,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        eq_( 0, ret )
        print out
        print open('outdir/tests.std.log').read()
        self._ensure_expected_output_files( 'outdir', 'tests' )

    def test_outdir_not_exist( self ):
        assert not isdir( 'outdir' )
        out,ret = self._run_runsample( self.reads_by_sample, self.ref, 'tests', 'outdir' )
        print out
        print open('outdir/tests.std.log').read()
        eq_( 0, ret )
        assert isdir( 'outdir' )
        self._ensure_expected_output_files( 'outdir', 'tests' )

    def test_missing_executables_exit_1( self ):
        # Change path so tempdir is first in lookup and then duplicate
        # some of the pipeline scripts and just have them return 1
        with open(join(self.tempdir,'samtools'),'w') as fh:
            fh.write( '#!/bin/bash\nexit 1\n' )
        os.chmod( join(self.tempdir,'samtools'), 0755 )
        import subprocess
        script = join( dirname( dirname( abspath( __file__ ) ) ), 'runsample.py' )
        cmd = 'export PATH={}:$PATH; {} {} {} {} -od {}'.format( self.tempdir, script, self.reads_by_sample, self.ref, 'tests', 'outdir' )
        ret = subprocess.call( cmd, shell=True )
        assert ret != 0, "Return code was 0 even though some executables returned 1"

    def test_ensure_proper_log( self ):
        # Just check that logfile gets stuff in it
        outdir = 'outdir'
        project = 'project'
        logfile = join( outdir, project+'.log' )
        out,ret = self._run_runsample( self.reads_by_sample, self.ref, project, outdir )
        loglines = None
        with open( logfile ) as fh:
            loglines = fh.read().splitlines()
        # 5 stages + start/finish
        ok_( 7 <= len(loglines), "Should be at least 7 loglines in log file" )
