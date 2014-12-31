import json

import tempdir

from imports import *
from test_bqd import Base as Base_

class Base(Base_):
    modulepath = 'ngs_mapper.coverage'

class TestFilterRefs(Base):
    functionname = 'filter_refs'

    def setUp(self):
        super(TestFilterRefs,self).setUp()
    
        # Tests
        self.refs = ['aa','b','caa','daa']

    def test_include_and_exclude(self):
        f = self._C(self.refs, ['aa'], ['aa'])
        eq_(f, set(['aa', 'caa','daa']))

    def test_no_filter(self):
        f = self._C(self.refs, [], ['unmapped_reads'])
        eq_(f, set(self.refs))

    def test_include_only(self):
        f = self._C(self.refs, ['aa'], [])
        eq_(f, set(['aa','caa','daa']))

    def test_exclude_only(self):
        f = self._C(self.refs, [], ['aa'])
        eq_(f, set(['b']))

@patch('__builtin__.open',Mock())
@patch('ngs_mapper.coverage.glob')
@patch('ngs_mapper.coverage.json')
class TestLoadProjectQualdepth(Base):
    functionname = 'load_project_qualdepth'

    def setUp(self):
        super(TestLoadProjectQualdepth,self).setUp()
        self.json = {'Ref1':1, 'Ref2':2}

    def test_loads_json(self, mock_json, mock_glob):
        mock_json.load.return_value = self.json
        mock_glob.return_value = ['']
        r = self._C('')
        ok_(isinstance(r, dict),'Did not return dict')

    @raises(ValueError)
    def test_project_missing_qualdepthfile(self, mock_json, mock_glob):
        mock_json.load.return_value = self.json
        mock_glob.return_value = []
        self._C('')

class TestRefsFromProject(Base):
    functionname = 'refs_from_project'

    @raises(ValueError)
    @patch('ngs_mapper.coverage.glob')
    def test_no_qualdepth_found(self, mock_glob):
        mock_glob.return_value = []
        self._C('', [], [])

    @patch('ngs_mapper.coverage.load_project_qualdepth')
    def test_returns_set(self, mock_lpqd):
        mock_lpqd.return_value = {'Ref1':1,'Ref2':2,'Ref1':3}
        r = self._C('', [], [])
        ok_(isinstance(r,set),'Did not return set')

class TestGetAllRefs(Base):
    functionname = 'get_allrefs'

    @patch('ngs_mapper.coverage.load_project_qualdepth')
    def test_returns_set(self, mock_lpqd):
        mock_lpqd.return_value = {'Ref1':1,'Ref2':2,'Ref1':3}
        r = self._C(['','',''], [], [])
        ok_(isinstance(r,set),'Did not return set')
        eq_(set(['Ref1','Ref2']),r)

class TestGetPerreferenceFromProjects(Base):
    functionname = 'get_perreference_from_projects'

    def setUp(self):
        super(TestGetPerreferenceFromProjects,self).setUp()
        from ngs_mapper.bqd import REGIONTYPES

        self.allrefs = ['Ref1','Ref2','Ref3']
        self.refax = {ref:Mock() for ref in self.allrefs}
        self.lineargs = self._make_lineargs(REGIONTYPES)
        self.qualdepth = {
            ref:{
                'depths':[0*5,10*5,10*5],
                'avgquals':[0*5,40*5,0*5],
                'length': 15,
                'minq': 0,
                'maxq': 40,
                'mind': 0,
                'maxd': 10,
            } for ref in self.allrefs
        }

    @raises(ValueError)
    @patch('ngs_mapper.coverage.glob',Mock(return_value=[]))
    def test_missing_qualdepth(self):
        self._C([''], [], [], 0, 0, 0, {})

    @patch('ngs_mapper.coverage.load_project_qualdepth')
    def test_returns_dict(self, mock_lpqd):
        mock_lpqd.return_value = self.qualdepth
        projects = ['1','2','3']
        r = self._C(projects, self.allrefs,
             self.refax, 0, 0, 0, self.lineargs
        )
        ok_(isinstance(r,dict),'Did not return dictionary')

class TestSetFigureSize(Base):
    functionname = 'set_figure_size'

    def setUp(self):
        super(TestSetFigureSize,self).setUp()

        self.perreference = {
            'Ref1': [[]],
            'Ref2': [[]],
        }

    def test_limits_size_to_400_inches(self):
        fig = Mock()
        fig.get_dpi.return_value = 72
        self.perreference['Ref1'] = [[1]]*10000
        self._C(self.perreference, fig)
        ok_(fig.set_dpi.called, 'Did not call set_dpi')
        ok_(fig.set_dpi.call_args_list[0][0][0] != 72)
        args = fig.set_size_inches.call_args
        ok_(args[0][1] == 400, 'Did not set length to 400')

    def test_sets_size_inches(self):
        fig = Mock()
        self._C(self.perreference, fig)
        ok_(fig.set_size_inches.called, 'Did not call set_size_inches')

    def test_numrefs_lt2_creates_graphic_correctly(self):
        from mock import call
        fig = Mock()
        del self.perreference['Ref2']
        self._C(self.perreference, fig)
        cl = fig.set_size_inches.call_args_list
        eq_([call(20.0,2)], cl)

class PlotBase(Base):
    def setUp(self):
        super(PlotBase,self).setUp()
        self.samplenames = ['1','2']
        self.line2ds = [
            self._mock_line2d([0,0],[i,i])
            for i in range(len(self.samplenames))
        ] + [self._mock_line2d([0,0],[0,0])]
        self.perreference = {
            'Ref1': [self.line2ds],
            'Ref2': [self.line2ds],
        }
        self.refax = {
            ref:Mock() for ref in self.perreference
        }
        self.projects = [p for p in self.samplenames]

    def _mock_line2d(self, xdata, ydata):
        class Line2D(object):
            _ydata = []
            _xdata = []
            def set_ydata(self,y):
                self._ydata = y
            def get_ydata(self,x):
                self._xdata = x
            def get_ydata(self):
                return self._ydata
            def get_xdata(self):
                return self._xdata
        line = Line2D()
        line._ydata = ydata
        line._xdata = xdata
        return line

class TestPlotReference(PlotBase):
    functionname = 'plot_reference'

    def test_adds_lines(self):
        ax = Mock()
        r = self._C('Ref1', self.perreference['Ref1'], self.samplenames, ax)
        eq_( [''] + self.samplenames, r )
        eq_(len(self.line2ds), ax.add_line.call_count)

class TestPlotAllReferences(PlotBase):
    functionname = 'plot_all_references'

    def test_plots_all_references(self):
        self._C(self.projects, self.perreference, self.refax)

class TestCreateLegend(Base):
    functionname = 'create_legend'

    @patch('ngs_mapper.coverage.Line2D')
    def test_correct_legend(self, mock_line2d):
        from ngs_mapper.bqd import REGIONTYPES
        fig = Mock()
        lineargs = self._make_lineargs(REGIONTYPES)
        mock_line2d.side_effect = range(len(REGIONTYPES))
        self._C(fig, REGIONTYPES, lineargs)
        ok_(fig.legend.called, 'Did not call figure.legend')
        
        eq_( [0,1,2,3,4], fig.legend.call_args_list[0][0][0] )

class TestCreateFigureForProjects(Base):
    functionname = 'create_figure_for_projects'

    @patch('ngs_mapper.coverage.plt')
    @patch('ngs_mapper.coverage.get_allrefs')
    @patch('ngs_mapper.coverage.get_perreference_from_projects')
    @patch('ngs_mapper.coverage.set_figure_size')
    @patch('ngs_mapper.coverage.plot_all_references')
    @patch('ngs_mapper.coverage.create_legend')
    @patch('ngs_mapper.coverage.gridspec')
    def test_returns_figure(self,*args):
        from ngs_mapper.bqd import REGIONTYPES
        r = self._C([''], [], [], self._make_lineargs(REGIONTYPES), [0,25,10])
        eq_( args[-1].figure(), r )

class TestMain(Base):
    functionname = 'main'
    
    @patch('ngs_mapper.coverage.gridspec',MagicMock())
    @patch('ngs_mapper.coverage.parse_args')
    @patch('ngs_mapper.coverage.plt')
    def test_runs(self, mock_plt, mock_args):
        args = Mock()
        mock_args.return_value = args

        args.lowcov = 10
        args.lowqual = 25
        args.gap = 0
        args.exclude = []
        args.include = []
        args.output = 'output.png'

        with tempdir.TempDir() as t:
            projdir = join(t,'sample')
            args.projects = [projdir]
            qdepthfile = join(projdir,'sample.bam.qualdepth.json')
            os.mkdir(projdir)
            ref1qd = self._make_qualdepth()
            ref2qd = self._make_qualdepth()
            qd = {'unmapped_reads':0, 'Ref1':ref1qd, 'Ref2':ref2qd}
            with open(qdepthfile,'w') as fh:
                json.dump(qd, fh)
            self._C()
