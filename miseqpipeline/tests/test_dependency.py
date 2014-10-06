from imports import *
import tempdir

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
            sys.stdout.write( 'Note: checking out \'{0}\'.\n'.format(
                    cmd[2]
                )
            )
            sys.stdout.write( 'HEAD is now at {0}... init\n'.format(
                    cmd[2]
                )
            )
        else:
            return 1
    elif cmd[0] == 'make':
        if isfile('Makefile'):
            make_mock_exec_file('bwa')
            make_mock_exec_file('samtools')
            make_mock_exec_file('bcftools/bcftools')
            sys.stdout.write('Making bwa\n')
            sys.stdout.write('Making samtools\n')
            sys.stdout.write('Making bcftools/bcftools\n')
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

    distros = [
        ('Red Hat Enterprise Linux Workstation', '6.5', 'Santiago'),
        ('CentOS', '6.5', 'Final'),
        ('Ubuntu', '14.04', 'Precise'),
        ('debian', 'wheezy/sid', ''),
    ]

    def setUp( self ):
        super(Base,self).setUp()
        self.prefix = 'prefixdir'
        self.bindir = join(self.prefix,'bin')
        self.libdir = join(self.prefix,'lib')
        self.bwapath = join(self.bindir,'bwa')
        self.sampath = join(self.bindir,'samtools')
        self.bcfpath = join(self.bindir,'bcftools')
        self.trimmopath = join(self.libdir, 'Trimmomatic-0.01', 'trimmomatic-0.01.jar')

    def _print_prefix( self, prefix ):
        for root, dirs, files in os.walk(prefix):
            for f in files:
                print join(root,f)

    def _exist( self, path ):
        ok_(exists(path), "{0} did not exist".format(path))

    def _exist_executable( self, path ):
        self._exist( path )
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
class TestVerifySamtoolsInstall(Base):
    functionname = 'verify_samtools_install'

    def setUp(self):
        super(TestVerifySamtoolsInstall,self).setUp()
        os.makedirs(self.bindir)

    def test_executable_exists_in_bin_not_executable_returns_false(self):
        open(self.sampath,'w').close()
        open(self.bcfpath,'w').close()
        os.chmod(self.sampath,0644)
        os.chmod(self.bcfpath,0644)
        
        r = self._C(self.prefix)
        eq_( False, r )

    def test_executable_exists_in_bin_executable_returns_true(self):
        open(self.sampath,'w').close()
        open(self.bcfpath,'w').close()
        os.chmod(self.sampath,0755)
        os.chmod(self.bcfpath,0755)
        
        r = self._C(self.prefix)
        eq_( True, r )

    def test_integration_verifies_install_samtools(self):
        from miseqpipeline.dependency import install_samtools
        install_samtools(self.samtools_github_url, 'HEAD', self.prefix)
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

class TestVerifyExistExecutable(Base):
    functionname = 'prefix_has_files'

    def setUp(self):
        super(TestVerifyExistExecutable,self).setUp()
        self.exefile = join(self.bindir,'myexe')
        self.libfile = join(self.libdir,'libfile')
        self.nested = join(self.libdir,'mylib','libfile')
        make_mock_exec_file( self.exefile )
        make_mock_exec_file( self.libfile )
        make_mock_exec_file( self.nested )
        self.files = [
            'bin/myexe',
            'lib/libfile',
            'lib/mylib/libfile',
        ]

    def test_verifies_paths_in_list_exist_and_executable(self):
        r = self._C( self.prefix, self.files )
        eq_( [], r )

    def test_missing_files_are_returned_as_list(self):
        os.unlink(self.exefile)
        r = self._C( self.prefix, self.files )
        eq_( self.files[0:1], r )
        os.unlink(self.libfile)
        r = self._C( self.prefix, self.files )
        eq_( self.files[0:2], r )
        os.unlink(self.nested)
        r = self._C( self.prefix, self.files )
        eq_( self.files, r )

    def test_checks_accessmode_for_files(self):
        modes = [os.X_OK, os.F_OK, os.F_OK]
        os.chmod(self.exefile, 0644)
        r = self._C( self.prefix, zip(self.files,modes) )
        eq_( r, ['bin/myexe'] )

def urllib_urlretrieve_mock( url, filename=None, reporthook=None, data=None ):
    ''' Just creates a compressed text file with hello in it '''
    if filename is None:
        dst = basename(url)
    else:
        dst = filename

    base_name, ext = splitext(basename(url))
    # Mocked trimmomatic zip
    os.mkdir('Trimmomatic-0.01')
    with open('Trimmomatic-0.01/trimmomatic-0.01.jar','w') as fh:
        fh.write('#!/bin/bash\n')
        fh.write('echo "trimmomatic"\n')
    with open('hello.txt','w') as fh:
        fh.write('hello\n')
    if ext == '.zip':
        import zipfile
        with zipfile.ZipFile(dst,'w') as zfh:
            zfh.write('Trimmomatic-0.01/trimmomatic-0.01.jar')
            zfh.write('hello.txt')
    else:
        import tarfile
        compression = ''
        if ext in ('.gz','.gzip','.tgz'):
            compression = 'w:gz'
        elif ext in ('.bz','.bzip','.bz2','.tbz'):
            compression = 'w:bz2'
        elif ext == '.tar':
            compression = 'w'
        else:
            raise Exception('Cannot determine compression for file {0}'.format(dst))
        with tarfile.open(dst, compression) as tar:
            tar.add('Trimmomatic-0.01')
            tar.add('hello.txt')
    shutil.rmtree('Trimmomatic-0.01')
    os.unlink('hello.txt')
    return (dst, '')

@patch('miseqpipeline.dependency.urllib',Mock(urlretrieve=urllib_urlretrieve_mock))
class TestDownloadUnpack(Base):
    functionname = 'download_unpack'

    def setUp(self):
        super(TestDownloadUnpack,self).setUp()
        self.unpackdir = 'unpackdir'
        self.hellopath = join(self.unpackdir,'hello.txt')

    def _check_hellofile( self, hellopath ):
        ok_( exists(hellopath), 'Did not unpack {0}'.format(hellopath) )
        with open(hellopath) as fh:
            eq_( 'hello\n', fh.read(), 'Contents of hello.txt are not correct' )

    def test_unpacks_zip_files(self):
        self._C( 'http://www.example.com/zipfile.zip', self.unpackdir )
        self._check_hellofile(self.hellopath)

    def test_unpacks_tar_files(self):
        self._C( 'http://www.example.com/zipfile.tar', self.unpackdir )
        self._check_hellofile(self.hellopath)

    def test_unpacks_tgz_files(self):
        self._C( 'http://www.example.com/zipfile.tgz', self.unpackdir )
        self._check_hellofile(self.hellopath)

    def test_unpacks_tbz_files(self):
        self._C( 'http://www.example.com/zipfile.tar.bz2', self.unpackdir )
        self._check_hellofile(self.hellopath)

    def test_unpacks_into_abspath(self):
        self._C( 'http://www.example.com/zipfile.tar.bz2', abspath(self.unpackdir) )
        self._check_hellofile(self.hellopath)

    def test_unpacks_into_relativepath(self):
        self._C( 'http://www.example.com/zipfile.tar.bz2', self.unpackdir )
        self._check_hellofile(self.hellopath)

    @raises(ValueError)
    def test_unkown_format_raises_exception(self):
        self._C( 'http://www.example.com/unkown.txt', '.' )

    def test_ensure_no_leftover_files_in_current_directory(self):
        self._C( 'http://www.example.com/file.zip', self.unpackdir)
        # These files are ok for our tests
        okfiles = [self.unpackdir]
        eq_( sorted(okfiles), sorted(os.listdir('.')) )

@patch('miseqpipeline.dependency.urllib',Mock(urlretrieve=urllib_urlretrieve_mock))
class TestInstallTrimmomatic(Base):
    functionname = 'install_trimmomatic'

    def test_creates_intermediate_directories(self):
        self._C( self.trimmomatic_download_url, self.libdir )
        self._exist( self.trimmopath )

    def test_intermediate_directories_exist(self):
        os.makedirs(self.libdir)
        self._C( self.trimmomatic_download_url, self.libdir )
        self._exist( self.trimmopath )

    def test_absolute_path_to_prefix(self):
        self.libdir = abspath(self.libdir)
        self._C( self.trimmomatic_download_url, self.libdir )
        self._exist( self.trimmopath )

class TestVerifyTrimmomatic(Base):
    functionname = 'verify_trimmomatic'

    def _make_trimmo_version( self, prefix, version ):
        ''' Make trimmodir and trimmo jar with version '''
        trimname = 'Trimmomatic-{0}'.format(version)
        trimodir = join(self.libdir,trimname)
        trimmopath = join(trimodir,trimname.lower()+'.jar')
        os.makedirs(trimodir)
        with open(trimmopath,'w') as fh:
            fh.write('i exist')
        return trimodir

    def test_jarfile_missing_from_dstprefix_lib_trimmomatic(self):
        os.makedirs(self.libdir)
        self._print_prefix(self.prefix)

        r = self._C( self.prefix )
        eq_( False, r )

    def test_jarfile_exists_in_dstprefix_lib_trimmomatic(self):
        os.makedirs(self.libdir)
        self._make_trimmo_version(self.prefix, '0.01')
        self._print_prefix(self.prefix)

        r = self._C( self.prefix, '0.01' )
        eq_( True, r )

    def test_incorrect_version_returns_false(self):
        os.makedirs(self.libdir)
        self._make_trimmo_version(self.prefix, '0.01')
        self._print_prefix(self.prefix)

        r = self._C(self.prefix, '0.32')
        eq_(False, r)

    def test_correct_version_with_multiple_installed_returns_true(self):
        os.makedirs(self.libdir)
        self._make_trimmo_version(self.prefix, '0.01')
        self._make_trimmo_version(self.prefix, '0.32')
        self._print_prefix(self.prefix)

        r = self._C(self.prefix, '0.32')
        eq_( True, r )

class TestGetDistributionPackageManager(Base):
    functionname = 'get_distribution_package_manager'

    def test_returns_correct_distribution_with_version(self):
        pkgmanagers = [
            'yum',
            'yum',
            'apt-get',
            'apt-get',
        ]
        for dist, pkgmanager in zip(self.distros,pkgmanagers):
            with patch('miseqpipeline.dependency.platform') as platform:
                platform.linux_distribution.return_value = dist
                r = self._C()
                eq_( pkgmanager, r )

    def test_raises_exception_with_unknown_distrubution(self):
        from miseqpipeline.dependency import UnknownDistributionError
        with patch('miseqpipeline.dependency.platform') as platform:
            platform.linux_distribution.return_value = (
                'Unknown Distribution', '1.0', 'Final'
            )
            try:
                r = self._C()
                ok_( False, "Did not raise UnknownDistributionError" )
            except UnknownDistributionError as e:
                ok_( True )

@patch('miseqpipeline.dependency.os')
@patch('miseqpipeline.dependency.subprocess')
@patch('miseqpipeline.dependency.platform')
class TestInstallSystemPackages(Base):
    functionname = 'install_system_packages'

    def test_installs_yum_packages(self, platform, subprocess, os):
        os.getuid.return_value = 0
        platform.linux_distribution.return_value = self.distros[1]
        subprocess.check_call.return_value = 0
        r = self._C( ['pkg1', 'pkg2'] )
        subprocess.check_call.assert_has_calls(
            call(['yum','install','-y','pkg1','pkg2'])
        )

    def test_installs_aptget_packages(self, platform, subprocess, os):
        os.getuid.return_value = 0
        platform.linux_distribution.return_value = self.distros[2]
        subprocess.check_call.return_value = 0
        r = self._C( ['pkg1', 'pkg2'] )
        subprocess.check_call.assert_has_calls(
            call(['apt-get','install','-y','pkg1','pkg2'])
        )

    def test_unknown_distribution_raises_unknowndistribution(self, platform, subprocess, os):
        from miseqpipeline.dependency import UnknownDistributionError
        os.getuid.return_value = 0
        platform.linux_distribution.return_value = ('Unknown','14.04','Precise')
        subprocess.check_call.side_effect = OSError
        try:
            r = self._C(['pkg1'])
            ok_( False, "Did not raise UnknownDistributionError" )
        except UnknownDistributionError as e:
            ok_(True)

    def test_uid_is_not_zero_raises_exception(self, platform, subprocess, os):
        from miseqpipeline.dependency import UserNotRootError
        os.getuid.return_value = 500
        platform.linux_distribution.return_value = self.distros[2]
        subprocess.check_call.side_effect = subprocess.CalledProcessError

        try:
            r = self._C( ['pkg1'] )
            ok_(False, "Did not raise Exception when user was not root")
        except UserNotRootError as e:
            ok_(True)

@attr('current')
@patch('miseqpipeline.dependency.platform')
class TestGetDistributionPackageList(Base):
    functionname = 'get_distribution_package_list'

    def setUp(self):
        super(TestGetDistributionPackageList, self).setUp()
        self.pkglist = {
            'yum': [
                'yum-pkg1',
                'yum-pkg2'
            ],
            'apt-get': [
                'aptget-pkg1',
                'aptget-pkg2'
            ]
        }
        self.pkglistfile = self._write_pkglistfile(
            self.pkglist,'system_packages.lst'
        )

    def _write_pkglistfile( self, pkglistdict, outpath ):
        import json
        with open(outpath,'w') as fh:
            json.dump(pkglistdict, fh, indent=4)
        return outpath

    def test_gets_correct_list_for_yum_platform(self, platform):
        platform.linux_distribution.return_value = self.distros[0]
        r = self._C( self.pkglistfile )
        eq_( self.pkglist['yum'], r )

    def test_gets_correct_list_for_aptget_platform(self, platform):
        platform.linux_distribution.return_value = self.distros[2]
        r = self._C( self.pkglistfile )
        eq_( self.pkglist['apt-get'], r )

    def test_pkglistfile_does_not_contain_pkgmanager_entry_raises_exception(self, platform):
        from miseqpipeline.dependency import MissingPackageListEntry
        platform.linux_distribution.return_value = self.distros[0]

        # Remove yum entry
        del self.pkglist['yum']
        self.pkglistfile = self._write_pkglistfile(
            self.pkglist, 'system_packages.lst'
        )

        try:
            r = self._C( self.pkglistfile )
            ok_(False, "Did not raise exception")
        except MissingPackageListEntry as e:
            ok_(True)
