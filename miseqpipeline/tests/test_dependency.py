from imports import *

def mock_bwa_subprocess_call( *args, **kwargs ):
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
        os.makedirs( join(clonedir,'.git') )
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
            with open('bwa','w') as fh:
                fh.write('#!/bin/bash\n')
                fh.write('echo "BWA"\n')
            os.chmod('bwa',0755)
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

@patch('miseqpipeline.dependency.subprocess',Mock(call=mock_bwa_subprocess_call))
class TestInstallBWA(Base):
    functionname = 'install_bwa'

    def test_runs_from_tempdir_and_removes_it(self):
        curlisting = os.listdir('.')
        tmplisting = os.listdir('/tmp')

        # Just run it 
        self._C(self.bwa_github_url, '', os.getcwd())

        # Current directory and temp should be the same
        # except if something wrote into temp in the mean time...fail
        eq_( curlisting + ['bin'], os.listdir('.') )
        eq_( tmplisting, os.listdir('/tmp') )

        ok_( os.access('bin/bwa', os.X_OK), 'bwa was not executable' )

    def test_installs_bwa_into_dstprefix(self):
        dstprefix = 'dstprefix'
        os.makedirs(dstprefix)

        self._C(self.bwa_github_url, '', dstprefix)

        expectedbwa = join(dstprefix,'bin','bwa')
        ok_(
            exists(expectedbwa),
            'Did not create bwa executable at {0}'.format(expectedbwa)
        )
        ok_(
            os.access(expectedbwa,os.X_OK),
            '{0} is not executable'.format(expectedbwa)
        )

    def test_bwa_executable_already_exists_and_is_overwritten(self):
        os.mkdir('bin')
        open('bin/bwa','w').close()
        size = os.stat('bin/bwa').st_size
        
        self._C(self.bwa_github_url,'',os.getcwd())

        newsize = os.stat('bin/bwa').st_size

        ok_( size != newsize, 'Did not replace existing bwa' )
