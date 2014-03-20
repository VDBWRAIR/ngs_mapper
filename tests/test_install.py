from imports import *
from . import tdir

class Base(common.BaseClass):
    vpath = join( tdir, '.miseqpipeline' )

    @classmethod
    def setUpClass( klass ):
        super(Base,klass).setUpClass()
        # Only install once because it takes a long time
        klass.returncode, klass.output = klass.run_installer()

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
        )
        self.python_packages = (
            'matplotlib',
            'vcf',
            'numpy',
            'bwa' 
        )
        os.chdir(tdir) 

    @classmethod
    def run_installer( klass ):
        # Install to test package tempdir
        return klass.run_script( 'install.sh', klass.vpath )

class TestFunctional( Base ):
    def test_install_ran_successfully( self ):
        eq_( 0, self.returncode )

    def test_links_scripts( self ):
        binpath = join( self.vpath, 'bin' )
        print os.listdir( binpath )
        for script in self.scripts:
            path = join( binpath, script )
            try:
                os.stat(path)
            except:
                ok_( False, "Script {} was not linked into virtenv bin directory to path {}".format(script,path) )

    def test_creates_virtenv( self ):
        print os.listdir( '.' )
        for me in self.must_exist:
            ok_( exists( me ), "install did not create {}".format(me) )

    def test_uninstall_works( self ):
        print self.run_script( 'uninstall.sh' )
        for me in self.must_exist:
            ok_( not exists( me ), "uninstall did not remove {}".format(me) )

    def test_python_packages_importable( self ):
        # Activate the env
        activate_this = join( self.vpath, 'bin', 'activate_this.py' )
        execfile(activate_this, dict(__file__=activate_this))
        for pkg in self.python_packages:
            ok_( __import__(pkg), "Could not import {}".format(pkg) )
