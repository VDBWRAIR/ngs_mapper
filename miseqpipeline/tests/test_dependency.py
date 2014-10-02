from imports import *

def make_mock_exec_file( path ):
    ''' just make a mock bwa to be called '''
    path = abspath(path)
    # Ensure directories all the way up to basename are made
    if not isdir(dirname(path)):
        os.makedirs(dirname(path))
    # Just make a simple bash script that is executable
    with open(path,'w') as fh:
        fh.write('#!/bin/bash\n')
        fh.write('exit 0\n')
    os.chmod(path,0755)

def mock_bwasamtools_subprocess_call( *args, **kwargs ):
    ''' Mock all subprocess.call calls for bwa install '''
    cmd = args[0]
    if ' '.join(cmd[0:2]) == 'git clone':
        # SHould get the last part of the url which git will
        # use as the local dir to clone into
        clonedir = basename(cmd[2])
        # Already exists
        if isdir(clonedir):
            sys.stderr.write('fatal: destination path \'{0}\' already ' \
                'exists and is not an empty directory.\n'.format(
                    clonedir
                )
            )
            return 128
        # make that clonedir
        os.makedirs(join(clonedir,'.git'))
        sys.stdout.write('Cloning into \'{0}\'...\n'.format(clonedir))
        sys.stdout.write('done.\n')
        return 0
    elif ' '.join(cmd[0:2]) == 'git checkout':
        # Have to be in a git directory
        if isdir('.git'):
            open('Makefile','w').close()
        else:
            return 1
    elif cmd[0] == 'make':
        if isfile('Makefile'):
            make_mock_exec_file('bwa')
            make_mock_exec_file('samtools')
            make_mock_exec_file('bcftools/bcftools')
        else:
            sys.stderr.write('make: *** No targets specified and no ' \
                'makefile found.  Stop.\n')
            return 2
    else:
        raise Exception('Umm what?')

class Base(common.BaseClass):
    modulepath = 'miseqpipeline.dependency'

    bwa_github_url = 'https://github.com/lh3/bwa'
    samtools_github_url = 'https://github.com/samtools/samtools'
    trimmomatic_download_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip'

    def setUp( self ):
        super(Base,self).setUp()
        self.prefix = 'prefixdir'
        self.bindir = join(self.prefix,'bin')
        self.bwapath = join(self.bindir,'bwa')
        self.sampath = join(self.bindir,'samtools')
        self.bcfpath = join(self.bindir,'bcftools')

    def _exist_executable( self, path ):
        ok_(exists(path), "{0} did not exist".format(path))
        ok_(os.access(path,os.X_OK), "{0} is not executable".format(path))

@patch('miseqpipeline.dependency.subprocess',Mock(call=mock_bwasamtools_subprocess_call))
class TestInstallBWA(Base):
    functionname = 'install_bwa'

    def test_dstprefix_bin_does_not_exist_and_is_created(self):
        self._C(self.bwa_github_url,'HEAD',self.prefix)
        self._exist_executable(self.bwapath)

        ok_( isdir(self.bindir), "{0} not created".format(self.bindir) )

    def test_dstprefix_bin_exists(self):
        self._C(self.bwa_github_url,'HEAD',self.prefix)
        self._exist_executable(self.bwapath)

    def test_installs_bwa_into_dstprefix(self):
        self._C(self.bwa_github_url, '', self.prefix)
        expectedbwa = join(self.prefix,'bin','bwa')
        self._exist_executable(expectedbwa)

    def test_bwa_executable_already_exists_and_is_overwritten(self):
        os.mkdir('bin')
        open('bin/bwa','w').close()
        size = os.stat('bin/bwa').st_size
        
        self._C(self.bwa_github_url,'',os.getcwd())
        self._exist_executable('bin/bwa')

        newsize = os.stat('bin/bwa').st_size
        ok_( size != newsize, 'Did not replace existing bwa' )

@patch('miseqpipeline.dependency.subprocess',Mock(call=mock_bwasamtools_subprocess_call))
class TestVerifyBwaInstall(Base):
    functionname = 'verify_bwa_install'

    def setUp(self):
        super(TestVerifyBwaInstall,self).setUp()
        os.makedirs(self.bindir)

    def test_executable_exists_in_bin_not_executable_returns_false(self):
        open(self.bwapath,'w').close()
        os.chmod(self.bwapath,0644)
        
        r = self._C(self.prefix)
        eq_( False, r )

    def test_executable_exists_in_bin_executable_returns_true(self):
        open(self.bwapath,'w').close()
        os.chmod(self.bwapath,0755)
        
        r = self._C(self.prefix)
        eq_( True, r )

    def test_integration_verifies_install_bwa(self):
        from miseqpipeline.dependency import install_bwa
        install_bwa(self.bwa_github_url, 'HEAD', self.prefix)
        r = self._C(self.prefix)
        eq_( True, r )

@patch('miseqpipeline.dependency.subprocess',Mock(call=mock_bwasamtools_subprocess_call))
class TestInstallSamtools(Base):
    functionname = 'install_samtools'

    def _C( self, *args, **kwargs ):
        super(TestInstallSamtools,self)._C(*args, **kwargs)
        self._exist_executable(self.sampath)
        self._exist_executable(self.bcfpath)

    def test_dstprefix_bin_does_not_exist_and_is_created(self):
        self._C(self.samtools_github_url, 'HEAD', self.prefix)

    def test_dstprefix_bin_exists(self):
        os.makedirs(self.bindir)
        self._C(self.samtools_github_url, 'HEAD', self.prefix)

    def test_executable_already_exists_and_is_overwritten(self):
        os.makedirs(self.bindir)
        open(self.sampath,'w').close()
        open(self.bcfpath,'w').close()
        self._C(self.samtools_github_url, 'HEAD', self.prefix)

@patch('miseqpipeline.dependency.subprocess',Mock(call=mock_bwasamtools_subprocess_call))
class TestCloneCheckoutMakeCopy(Base):
    functionname = 'clone_checkout_make_copy'

    def _C( self, *args, **kwargs ):
        super(TestCloneCheckoutMakeCopy,self)._C(*args, **kwargs)
        for path in kwargs['copypaths']:
            execpath = join(self.bindir,basename(path))
            self._exist_executable(execpath)

    def test_ensures_prefix_bin_exists(self):
        copypaths = ['bwa','samtools','bcftools/bcftools']
        self._C('/path/to/myapp', 'HEAD', self.prefix, copypaths=copypaths)

    def test_prefix_bin_exists_no_error(self):
        os.makedirs(self.bindir)
        copypaths = ['bwa','samtools','bcftools/bcftools']
        self._C('/path/to/myapp', 'HEAD', self.prefix, copypaths=copypaths)

    def test_overwrites_existing_executable(self):
        os.makedirs(self.bindir)
        open(self.bwapath,'w').close()
        copypaths = ['bwa','samtools','bcftools/bcftools']
        self._C('/path/to/myapp', 'HEAD', self.prefix, copypaths=copypaths)

        content = open(self.bwapath).read()
        ok_(content != '', 'Did not overwrite existing executable')
