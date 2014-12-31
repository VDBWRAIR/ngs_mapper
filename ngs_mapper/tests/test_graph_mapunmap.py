from imports import *

from subprocess import check_output, CalledProcessError, STDOUT
import numpy as np

class Base(common.BaseBamRef):
    modulepath = 'ngs_mapper.graph_mapunmap'

class TestUnitArgs(Base):
    functionname = 'parse_args'
    script = 'graph_mapunmap.py'

    def setUp( self ):
        super(TestUnitArgs,self).setUp()

    def _mock_args( self ):
        args = Mock()
        args.outfile = ''
        args.jsons = []
        return args

    def test_outfile_unset( self ):
        res = self._C( [self.script, 'qd.json'] )
        # Should just be in current directory with default name
        eq_( 'mapunmap.png', res.outfile )

    def test_outfile_set( self ):
        res = self._C( [self.script, 'qd.json', '-o', 'mu.png'] )
        eq_( 'mu.png', res.outfile )

    def test_input_single( self ):
        res = self._C( [self.script, 'qd.json'] )
        eq_( ['qd.json'], res.jsons )

    def test_input_multiple( self ):
        res = self._C( [self.script, 'qd1.json', 'qd2.json'] )
        eq_( ['qd1.json','qd2.json'], res.jsons )

class TestUnitSampleFromFilename(Base):
    functionname = 'sample_from_filename'

    def test_returns_everything_before_bam( self ):
        res = self._C( 'sample_sample78-sample.bam' )
        eq_( 'sample_sample78-sample', res )

    @raises(ValueError)
    def test_filename_missing_bam( self ):
        res = self._C( 'samplename.json' )

    def test_filename_relabs( self ):
        res = self._C( '../test/samplename.bam.qualdepth.json' )
        eq_( 'samplename', res )

class TestUnitGetMapUnmap(Base):
    functionname = 'get_mapunmap'

    @patch('__builtin__.open',Mock())
    @patch('ngs_mapper.graph_mapunmap.sample_from_filename',Mock(side_effect=['S1','S2']))
    @patch('json.load')
    def test_gets_values( self, load ):
        load.side_effect = [
            {
                'unmapped_reads': 10,
                'Ref1': {
                    'mapped_reads': 100
                }
            },
            {
                'unmapped_reads': 20,
                'Ref1': {
                    'mapped_reads': 100
                },
                'Ref2': {
                    'mapped_reads': 100
                }
            }
        ]

        res = self._C( ['',''] )
        eq_( ['S1','S2'], res[0] )
        eq_( [100, 200], res[1] )
        eq_( [10, 20], res[2] )

class TestFunctional(Base):
    # Should make files with these extensions
    outfiles = ( '.mapunmap.png' )
    json = None

    @classmethod
    def setUpClass( klass ):
        super(TestFunctional,klass).setUpClass()
        from ngs_mapper.graphsample import make_json
        klass.json = make_json( klass.bam, klass.bam )

    def test_ensure_json( self ):
        # Just make sure the setUpClass is working
        import json
        assert json.load( open(self.json) )

    def _run_cmd( self, jsons, outpath=None ):
        script_path = 'graph_mapunmap.py'
        args = ' '.join( jsons )
        if outpath is not None:
            args += ' -o {}'.format(outpath)
        cmd = script_path + ' {} '.format(args)
        print "Running: {}".format(cmd)
        cmd = shlex.split( cmd )
        try:
            sout = check_output( cmd, stderr=STDOUT )
        except CalledProcessError as e:
            print e.output
            assert False
        return sout

    def test_creates_graphic_from_single_json( self ):
        res = self._run_cmd( [self.json] )
        assert isfile( 'mapunmap.png' )
        ok_( os.stat('mapunmap.png').st_size >= 30000 )

    def test_creates_graphic_from_multiple_json( self ):
        res = self._run_cmd( [self.json, self.json] )
        assert isfile( 'mapunmap.png' )
        ok_( os.stat('mapunmap.png').st_size >= 30000 )

    def test_outfile_works( self ):
        ofile = join( self.tempdir, 'out.png' )
        res = self._run_cmd( [self.json], ofile )
        assert isfile( ofile )
        ok_( os.stat(ofile).st_size >= 30000 )
