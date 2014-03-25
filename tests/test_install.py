from imports import *
from . import tdir

class Base(common.BaseClass):
    vpath = join( tdir, '.miseqpipeline' )

    @classmethod
    def setUpClass( klass ):
        super(Base,klass).setUpClass()
        # Uninstall first since __init__.py already installed for us
        # Only install once because it takes a long time
        klass.returncode, klass.output = klass.run_installer( klass.vpath )

    def setUp( self ):
        super( Base, self ).setUp()
        self.must_exist = (
            self.vpath,
            join(self.vpath,'bin'),
            join(self.vpath,'man1'),
            join(self.vpath,'lib','python2.7','site-packages')
        )
        # Must exist in vpath/bin
        self.scripts = (
            'base_caller.py',
            'consensuses.sh',
            'gen_flagstats.sh',
            'graph_mapunmap.py',
            'graphsample.py',
            'graphs.sh',
            'graph_times.py',
            'install.sh',
            'run_bwa_on_samplename.py',
            'runsample.py',
            'runsamplesheet.sh',
            'stats_at_refpos.py',
            'tagreads.py',
            'uninstall.sh',
            'vcf_consensus.py',
            'vcf_diff.py',
            'cutadapt'
        )
        self.python_packages = (
            'matplotlib',
            'vcf',
            'numpy',
            'bwa' 
        )
        os.chdir(tdir) 

    @classmethod
    def run_installer( klass, installpath ):
        # Install to test package tempdir
        script = klass.script_path('install.sh')
        return klass.run_script( '{} {}'.format( script, installpath ) )

class TestFunctional( Base ):
    def test_install_ran_successfully( self ):
        print self.output
        eq_( 0, self.returncode )
        ok_( 'failed' not in self.output )

    def test_links_scripts( self ):
        binpath = join( self.vpath, 'bin' )
        print os.listdir( binpath )
        for script in self.scripts:
            path = join( binpath, script )
            try:
                os.stat(path)
            except:
                print self.output
                ok_( False, "Script {} was not linked into virtenv bin directory to path {}".format(script,path) )

    def test_creates_virtenv( self ):
        print os.listdir( '.' )
        for me in self.must_exist:
            ok_( exists( me ), "install did not create {}".format(me) )

    @timed(1)
    def test_wrong_python_version( self ):
        # Make a mock python that returns an older version number
        with open( 'python', 'w' ) as fh:
            fh.write( "#!/bin/bash\n" )
            fh.write( "echo 2.6.0" )
        # Make it executable
        st = os.stat('python')
        import stat
        os.chmod( 'python', st.st_mode | stat.S_IEXEC )
        script = self.script_path('install.sh')
        path=tdir + ':' + os.environ['PATH']
        ret, output = self.run_script( 'export PATH={}; {} {}'.format(path,script,self.vpath) )
        eq_( 'Please ensure you have python 2.7.3 or greater but less than 3', output.rstrip() )
        eq_( 1, ret )

    def test_uninstall_works( self ):
        print self.run_script( self.script_path( 'uninstall.sh' ) )
        for me in self.must_exist:
            ok_( not exists( me ), "uninstall did not remove {}".format(me) )

    def test_python_packages_importable( self ):
        # Activate the env
        activate_this = join( self.vpath, 'bin', 'activate_this.py' )
        execfile(activate_this, dict(__file__=activate_this))
        for pkg in self.python_packages:
            ok_( __import__(pkg), "Could not import {}".format(pkg) )
