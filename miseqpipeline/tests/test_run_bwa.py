from imports import *

class Base(BaseClass):
    modulepath = 'miseqpipeline.run_bwa'

class TestUnitBWAMem(Base):
    functionname = 'bwa_mem'

    def test_bwa_mem_nonpaired(self):
        with patch('miseqpipeline.run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=0 ) ) ) as b:
            with patch('miseqpipeline.run_bwa.index_ref',Mock(return_value=True)) as a:
                result = self._C( 'F.fq', mate=None, ref='ref.fna' )
                eq_( 'bwa.sai', result )

    def test_bwa_mem_paired(self):
        with patch('miseqpipeline.run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=0 ) ) ) as b:
            with patch('miseqpipeline.run_bwa.index_ref',Mock(return_value=True)) as a:
                result = self._C( 'F.fq', mate='R.fq', ref='ref.fna' )
                eq_( 'bwa.sai', result )

    def test_bwa_mem_output_arg(self):
        with patch('miseqpipeline.run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=0 ) ) ) as b:
            with patch('miseqpipeline.run_bwa.index_ref',Mock(return_value=True)) as a:
                result = self._C( 'F.fq', mate='R.fq', ref='ref.fna', output='file.sai' )
                eq_( 'file.sai', result )

    def test_bwa_mem_fails(self):
        with patch('miseqpipeline.run_bwa.BWAMem', return_value=Mock( run=Mock( return_value=1 ) ) ) as b:
            with patch('miseqpipeline.run_bwa.index_ref',Mock(return_value=True)) as b:
                result = self._C( 'F.fq', mate='R.fq', ref='ref.fna', output='file.sai' )
                eq_( 1, result )

    @patch('miseqpipeline.run_bwa.index_ref', Mock(return_value=False))
    def test_ref_index_fails(self):
        from miseqpipeline.run_bwa import bwa_mem, InvalidReference
        try:
            self._C( 'F.fq', mate='R.fq', ref='ref.fna', output='file.sai' )
        except InvalidReference as e:
            pass
        else:
            assert False, "Did not raise InvalidReference"

    @patch('miseqpipeline.run_bwa.index_ref')
    @patch('miseqpipeline.run_bwa.compile_refs')
    @patch('miseqpipeline.run_bwa.BWAMem')
    def test_ref_index_directory(self, bwamem_mock, compile_refs_mock, index_ref_mock ):
        bwamem_mock.return_value.run.return_value = 1
        compile_refs_mock.return_value = 'ref_compiled.fna'
        index_ref_mock.return_value = 1
        from miseqpipeline.run_bwa import InvalidReference
        ret = self._C( 'F.fq', mate='R.fq', ref='ref_compiled.fna', output='file.sai' )
        eq_( 1, ret )

    @patch('miseqpipeline.run_bwa.index_ref')
    @patch('miseqpipeline.run_bwa.BWAMem')
    @patch('miseqpipeline.run_bwa.which_bwa',Mock(return_value='bwa'))
    def test_set_threads(self, bwamem_mock,index_ref_mock ):
        bwamem_mock.return_value.run.return_value = 1
        index_ref_mock.return_value = 1
        from miseqpipeline.run_bwa import InvalidReference

        ret = self._C( 'F.fq', mate='R.fq', ref='ref.fna', output='file.sai', t=8 )
        bwamem_mock.assert_called_with( 'ref.fna', 'F.fq', 'R.fq', bwa_path='bwa', t=8 )

        ret = self._C( 'F.fq', ref='ref.fna', output='file.sai', t=8 )
        bwamem_mock.assert_called_with( 'ref.fna', 'F.fq', bwa_path='bwa', t=8 )

class TestUnitParseArgs(Base):
    functionname = 'parse_args'

    @raises(SystemExit)
    def test_ref_reads_required( self ):
        res = self._C( [] )

    def test_ref_reads_set( self ):
        res = self._C( ['fake_read', 'fake_ref'] )
        eq_( res.reads, 'fake_read' )
        eq_( res.reference, 'fake_ref' )

    def test_platform_select_none( self ):
        res = self._C( ['fake_read', 'fake_ref'] )
        eq_( res.platforms, ['MiSeq','Sanger','Roche454','IonTorrent'] )

    def test_platform_select_single( self ):
        res = self._C( ['fake_read', 'fake_ref', '--platforms', 'Sanger'] )
        eq_( res.platforms, ['Sanger'] )

    def test_output_path( self ):
        res = self._C( ['fake_read', 'fake_ref', '-o', 'out.bam'] )
        eq_( res.output, 'out.bam' )

    def test_output_path_long( self ):
        res = self._C( ['fake_read', 'fake_ref', '--output', 'out.bam'] )
        eq_( res.output, 'out.bam' )

    def test_output_path_default( self ):
        res = self._C( ['fake_read', 'fake_ref'] )
        eq_( res.output, 'bwa_mem.bam' )

    @raises(SystemExit)
    def test_invalid_platform( self ):
        res = self._C( ['fake_read', 'fake_ref', '--platforms', 'invalid'] )

    def test_keep_temp_defaultoff( self ):
        res = self._C( ['fake_read', 'fake_ref'] )
        eq_( res.keep_temp, False )

    def test_keep_temp_set( self ):
        res = self._C( ['fake_read', 'fake_ref', '--keep-temp'] )
        eq_( res.keep_temp, True )

    def test_threads_set( self ):
        res = self._C( ['fake_read', 'fake_ref', '-t', '5'] )
        eq_( res.threads, 5 )

# Pretty sure this isn't the way to do this, but I'm learning here
@patch('shutil.move')
@patch('shutil.rmtree')
@patch('miseqpipeline.run_bwa.parse_args')
@patch('miseqpipeline.bam.mergebams')
@patch('miseqpipeline.bam.indexbam')
@patch('miseqpipeline.bam.samtobam')
@patch('miseqpipeline.bam.sortbam')
@patch('miseqpipeline.run_bwa.bwa_mem')
@patch('miseqpipeline.run_bwa.compile_reads')
@patch('miseqpipeline.run_bwa.reads_by_plat')
@patch('miseqpipeline.run_bwa.compile_refs')
@patch('tempfile.mkdtemp')
class TestUnitMain(Base):
    functionname = 'main'

    def test_paired_readfiles(self,tmp_mock,ref_mock,reads_mock,compile_reads_mock, bwa_mem_mock, sort, convert, index, merge, parse_args, shrmtree, shmove):
        merge.side_effect = AssertionError("Should not merge single files")
        tmp_mock.return_value = 'tdir'
        os.mkdir('tdir')
        ref_mock.return_value = 'reference.fa'
        reads_mock.return_value = {'MiSeq':[('r1.fq','r2.fq')]}
        compile_reads_mock.return_value = {'F':'F.fq','R':'R.fq','NP':None}
        parse_args.return_value = Mock(reads='/reads', reference='/reference.fa', platforms=['MiSeq','Sanger'], keep_temp=False, threads=1)
        bwa_mem_mock.return_value = 'bwa_mem.bam'
        res = self._C()
        eq_( [call('F.fq','R.fq','/reference.fa','tdir/paired.sai',t=1)], bwa_mem_mock.call_args_list )
        eq_( [call([('r1.fq','r2.fq')],'tdir/reads')], compile_reads_mock.call_args_list )
        eq_( [call('/reads')], reads_mock.call_args_list )
        eq_( sort.call_count, 1 )
        eq_( convert.call_count, 1 )
        eq_( index.call_count, 1 )
        eq_( 1, shrmtree.call_count )

    def test_nonpaired_readfiles(self,tmp_mock,ref_mock,reads_mock,compile_reads_mock, bwa_mem_mock, sort, convert, index, merge, parse_args, shrmtree, shmove):
        merge.side_effect = AssertionError("Should not merge single files")
        tmp_mock.return_value = 'tdir'
        os.mkdir('tdir')
        ref_mock.return_value = 'reference.fa'
        reads_mock.return_value = {'Sanger':['r1.fq']}
        compile_reads_mock.return_value = {'F':None,'R':None,'NP':'NP.fq'}
        parse_args.return_value = Mock(reads='/reads', reference='/reference.fa', platforms=['MiSeq','Sanger'], keep_temp=False, threads=1)
        bwa_mem_mock.return_value = 'bwa_mem.bam'
        res = self._C()
        eq_( [call('NP.fq',ref='/reference.fa',output='tdir/nonpaired.sai',t=1)], bwa_mem_mock.call_args_list )
        eq_( [call(['r1.fq'],'tdir/reads')], compile_reads_mock.call_args_list )
        eq_( [call('/reads')], reads_mock.call_args_list )
        eq_( sort.call_count, 1 )
        eq_( convert.call_count, 1 )
        eq_( index.call_count, 1 )
        eq_( 1, shrmtree.call_count )
    
    def test_paired_and_nonpaired_readfiles(self,tmp_mock,ref_mock,reads_mock,compile_reads_mock, bwa_mem_mock, sort, convert, index, merge, parse_args, shrmtree, shmove):
        tmp_mock.return_value = 'tdir'
        os.mkdir('tdir')
        ref_mock.return_value = 'reference.fa'
        reads_mock.return_value = {'MiSeq':[('r1.fq','r2.fq')],'Sanger':['r3.fq']}
        compile_reads_mock.return_value = {'F':'F.fq','R':'R.fq','NP':'NP.fq'}
        parse_args.return_value = Mock(reads='/reads', reference='/reference.fa', platforms=['MiSeq','Sanger'], keep_temp=False, threads=1)
        bwa_mem_mock.return_value = 'bwa_mem.bam'
        res = self._C()
        eq_( [call('F.fq','R.fq','/reference.fa','tdir/paired.sai',t=1),call('NP.fq',ref='/reference.fa',output='tdir/nonpaired.sai',t=1)], bwa_mem_mock.call_args_list )
        eq_( [call([('r1.fq','r2.fq'),'r3.fq'],'tdir/reads')], compile_reads_mock.call_args_list )
        eq_( [call('/reads')], reads_mock.call_args_list )
        eq_( 2, sort.call_count )
        eq_( 2, convert.call_count )
        eq_( 1, index.call_count )
        eq_( 1, merge.call_count )
        eq_( 1, shrmtree.call_count )

    def test_keeptemp(self,tmp_mock,ref_mock,reads_mock,compile_reads_mock, bwa_mem_mock, sort, convert, index, merge, parse_args, shrmtree, shmove):
        tmp_mock.return_value = 'tdir'
        os.mkdir('tdir')
        shrmtree.side_effect = AssertionError("Should not remove files with keeptemp option")
        parse_args.return_value = Mock(reads='/reads', reference='/reference.fa', platforms=['MiSeq','Sanger'], keep_temp=True, threads=1)
        res = self._C()
        eq_( 0, shrmtree.call_count )

    def test_keeptemp_false(self,tmp_mock,ref_mock,reads_mock,compile_reads_mock, bwa_mem_mock, sort, convert, index, merge, parse_args, shrmtree, shmove):
        tmp_mock.return_value = 'tdir'
        os.mkdir('tdir')
        parse_args.return_value = Mock(reads='/reads', reference='/reference.fa', platforms=['MiSeq','Sanger'], keep_temp=False, threads=1)
        res = self._C()
        eq_( [call('tdir')], shrmtree.call_args_list )

    def test_utilizes_thread_arg(self,tmp_mock,ref_mock,reads_mock,compile_reads_mock, bwa_mem_mock, sort, convert, index, merge, parse_args, shrmtree, shmove):
        tmp_mock.return_value = 'tdir'
        os.mkdir('tdir')
        ref_mock.return_value = 'reference.fa'
        reads_mock.return_value = {'MiSeq':[('r1.fq','r2.fq')],'Sanger':['r3.fq']}
        compile_reads_mock.return_value = {'F':'F.fq','R':'R.fq','NP':'NP.fq'}
        bwa_mem_mock.return_value = 'bwa_mem.bam'
        parse_args.return_value = Mock(reads='reads', reference='reference.fa', platforms=['MiSeq','Sanger'], keep_temp=False, threads=8)
        res = self._C()
        bwa_mem_mock.assert_called_with('NP.fq', ref='reference.fa', output='tdir/nonpaired.sai', t=8)

    @attr('current')
    def test_bwa_error_should_raise_exception(self,*args):
        with patch('miseqpipeline.run_bwa.os') as os:
            from miseqpipeline.run_bwa import BWAError
            bwa_mem_mock = args[4]
            bwa_mem_mock.return_value = 1
            try:
                self._C()
                ok_(False,"Did not raise Exception for bwa error")
            except BWAError as e:
                ok_(True)

class TestIntegrateMainArgs(Base):
    def setUp(self):
        super(TestIntegrateMainArgs,self).setUp()
        os.mkdir( 'expected' )
        # Gives us easy access to some files to map against a reference and then
        # also have the expected bam files for paired, nonpaired and merged 
        # Keys of interest:
        #  NP - nonpaired
        #  P - paired
        #  REF - reference
        #  paired.bam
        #  nonpaired.bam
        #  merged.bam
        #  merged.bam.bai
        self.fixture_files = fixtures.unpack_integrated( 'expected' )

    def _CM( self, arglist ):
        # Workaround for the whole unittest using sys.argv
        from miseqpipeline.run_bwa import parse_args
        ns = parse_args( arglist )
        # Just ensure a static path for mkdtemp so we can ensure it exists or doesn't exist
        with patch( 'tempfile.mkdtemp' ) as mkdtemp:
            with patch( 'miseqpipeline.run_bwa.parse_args' ) as parse_args:
                mkdtemp.return_value = join(self.tempdir, 'tmpdir1')
                os.mkdir( join(self.tempdir,'tmpdir1') )
                parse_args.return_value = ns
                from miseqpipeline.run_bwa import main
                return main()

    def _eqsize( self, f1, f2, threshold=10000 ):
        # Just make sure file sizes are relatively similar
        # threshold size in bytes that they can differ
        s1 = os.stat( f1 ).st_size
        s2 = os.stat( f2 ).st_size
        diff = abs( s1 - s2 )
        assert diff <= threshold, "Size of {}({}) - {}({}) was greater than threshold of {}".format(
            f1, s1, f2, s2, threshold
        )

    def test_handles_sff( self ):
        sff = glob( join( fixtures.THIS, 'fixtures', 'reads', '*.sff' ) )[0]
        shutil.copy( sff, 'expected/reads' )
        sff = join( 'expected', 'reads', basename( sff ) )

        print "All files in expected directory:" + str( glob( 'expected/*' ) )
        ff = self.fixture_files
        argv = ['expected/reads', ff['REF']]
        r = self._CM( argv )
        assert os.stat( r )
        import subprocess
        out = subprocess.check_output( ['samtools', 'view', '{}'.format(r)] )
        rochecount = out.count( 'IA52U1' )
        print rochecount
        eq_( 100, rochecount, 'Sff file reads did not make it into bam file' )

        os.unlink( sff )

    def test_paired_and_nonpaired_get_merged(self):
        ff = self.fixture_files
        argv = ['expected/reads', ff['REF'], '--keep-temp']
        res = self._CM( argv )
        eq_( 'bwa_mem.bam', res )
        self._eqsize( ff['merged.bam'], res )
        self._eqsize( ff['merged.bam.bai'], res+'.bai' )
        assert os.path.exists( 'tmpdir1' ), "Did not keep temp directory"

    def test_paired_only(self):
        ff = self.fixture_files
        argv = ['expected/reads', self.fixture_files['REF'], '--platforms', 'MiSeq']
        res = self._CM( argv )
        eq_( 'bwa_mem.bam', res )
        self._eqsize( ff['paired.bam'], res )
        assert not os.path.exists( 'tmpdir1' ), "Temp directory still exists"
        assert os.path.exists( 'bwa_mem.bam.bai' )

    def test_nonpaired_only(self):
        ff = self.fixture_files
        print ff
        argv = ['expected/reads', self.fixture_files['REF'], '--platforms', 'Sanger']
        res = self._CM( argv )
        eq_( 'bwa_mem.bam', res )
        self._eqsize( ff['nonpaired.bam'], res )
        assert not os.path.exists( 'tmpdir1' ), "Temp directory still exists"
        assert os.path.exists( 'bwa_mem.bam.bai' )

    def test_output_path(self):
        ff = self.fixture_files
        argv = ['expected/reads', self.fixture_files['REF'], '-o', 'merged.bam']
        res = self._CM( argv )
        eq_( 'merged.bam', res )
        self._eqsize( ff['merged.bam'], res )
        assert not os.path.exists( 'tmpdir1' ), "Temp directory still exists"
        assert os.path.exists( 'merged.bam.bai' )

    def test_keepfiles(self):
        ff = self.fixture_files
        argv = ['expected/reads', self.fixture_files['REF'], '-o', 'merged.bam', '--keep-temp']
        res = self._CM( argv )
        eq_( 'merged.bam', res )
        self._eqsize( ff['merged.bam'], res )
        assert os.path.exists( 'tmpdir1' ), "Temp directory still exists"
        assert os.path.exists( 'merged.bam.bai' )

class TestFunctionalRunBWA(Base):
    functionname = 'bwa_mem'

    def setUp(self):
        super(TestFunctionalRunBWA,self).setUp()
        self.read1,self.read2,self.ref = fixtures.get_sample_paired_reads()
        self.sff = glob( join( fixtures.THIS, 'fixtures', 'reads', '*.sff' ) )[0]

    def test_maps_reads_paired(self):
        eq_( 'bwa.sai', self._C( self.read1, self.read2, ref=self.ref ) )
        assert exists( 'bwa.sai' ), "Did not create a sai file"
        assert os.stat('bwa.sai' ).st_size != 0, "sai file created is zero bytes"

    def test_maps_reads_single(self):
        eq_( 'bwa.sai', self._C( self.read1, ref=self.ref ) )
        assert exists( 'bwa.sai' ), "Did not create a sai file"
        assert os.stat('bwa.sai').st_size != 0, "sai file created is zero bytes"

    def test_output_param(self):
        eq_( 'file.sai', self._C( self.read1, ref=self.ref, output='file.sai' ) )
        assert exists( 'file.sai' ), "Did not create a sai file"
        assert os.stat('file.sai').st_size != 0, "sai file created is zero bytes"

    def test_ref_is_directory(self):
        import shutil
        os.mkdir( 'refs' )
        r1 = join('refs','ref1.fna')
        r2 = join('refs','ref2.fna')
        shutil.copy( self.ref, r1 )
        shutil.copy( self.ref, r2 )
        tot_size = os.stat(r1).st_size + os.stat(r2).st_size
        eq_( 'bwa.sai', self._C( self.read1, self.read2, ref='refs' ) ) 
        # bwa.bwa.compile_refs produces reference.fa inside of current directory
        eq_( tot_size, os.stat( 'reference.fa' ).st_size )
