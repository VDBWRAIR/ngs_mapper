from imports import *
from mock import mock_open, call
import yaml
from StringIO import StringIO

from ngs_mapper import config

class Base(common.BaseClass):
    modulepath = 'ngs_mapper.config'

    def setUp(self):
        super(Base,self).setUp()

        self.config = {
            'NGSDATA': self.tempdir,
            'base_caller': {
                'bias': {
                    'default': 10
                }
            }
        }

    def _create_yaml_from_config(self, config):
        import yaml
        return yaml.dump(config)

@patch('ngs_mapper.config.yaml', autospec=True)
class TestLoadConfigFile(Base):
    functionname = 'load_config'

    def test_loads_config_stream(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        config_stream = mock_open(read_data=self._create_yaml_from_config(self.config))()
        r = self._C(config_stream)
        eq_(self.tempdir, r['NGSDATA'])

    def test_loads_config_filepath(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        with patch('__builtin__.open') as mock_open:
            r = self._C('/path/to/file.yaml')
            eq_(self.tempdir, r['NGSDATA'])
            mock_open.assert_called_once_with('/path/to/file.yaml')
            mock_yaml.load.assert_called_once_with(mock_open())

    @patch('__builtin__.open', Mock())
    def test_invalid_config_raises_exception(self, mock_yaml):
        self.config['NGSDATA'] = '/missing/path'
        mock_yaml.load.return_value = self.config
        from ngs_mapper.config import InvalidConfigError
        assert_raises(InvalidConfigError, self._C, '/path/to/file.yaml')

    @attr('current')
    @patch('__builtin__.open', Mock())
    def test_returns_config_class_instance(self, mock_yaml):
        from ngs_mapper.config import Config
        mock_yaml.load.return_value = self.config
        r = self._C('/path/to/my.yaml')
        assert_is_instance(r, Config)

@attr('current')
class TestConfigClass(Base):
    def setUp(self):
        super(TestConfigClass,self).setUp()
        self.config = {
            'foo': 'bar'
        }
        self.inst = config.Config(self.config)

    def test_returns_values_from_getitem(self):
        eq_('bar', self.inst['foo'])

    def test_returns_values_from_attributes(self):
        eq_('bar', self.inst.foo)

    def test_raises_exception_missing_value_attribute(self):
        self.inst.yaml = {}
        assert_raises(config.InvalidConfigError, self.inst.__getattr__, 'foo')

    def test_raises_exception_missing_value_getitem(self):
        self.inst.yaml = {}
        assert_raises(config.InvalidConfigError, self.inst.__getitem__, 'foo')

class TestVerifyConfig(Base):
    functionname = 'verify_config'

    def test_bias_value_lt_1_raises_error(self):
        from ngs_mapper.config import InvalidConfigError
        self.config['base_caller']['bias']['default'] = 0
        assert_raises(InvalidConfigError, self._C, self.config)

    def test_ngsdata_not_set_raises_error(self):
        from ngs_mapper.config import InvalidConfigError
        self.config['NGSDATA'] = '/non/existant/path'
        assert_raises(InvalidConfigError, self._C, self.config)

    def test_valid_config_no_error(self):
        self.config['NGSDATA'] = self.tempdir
        self._C(self.config)

@patch('__builtin__.open', Mock())
@patch('pkg_resources.resource_stream')
@patch('ngs_mapper.config.yaml')
class TestLoadDefaultConfig(Base):
    functionname = 'load_default_config'

    def test_loads_from_pkg_resources(self, mock_yaml, mock_stream):
        mock_stream.return_value = mock_open(read_data=self._create_yaml_from_config(self.config))()
        mock_yaml.load.return_value = self.config
        r = self._C()
        eq_(self.tempdir, r['NGSDATA'])

@patch('ngs_mapper.config.yaml')
@patch('__builtin__.open', MagicMock())
class TestMakeExampleConfig(Base):
    functionname = 'make_example_config'

    def test_invalid_savepath_raises_exception(self, mock_yaml):
        assert_raises(ValueError, self._C, '/non/existant/path.yaml')

    def test_makes_config_cwd_is_default(self, mock_yaml):
        from ngs_mapper.config import load_config
        mock_yaml.load.return_value = self.config
        r = self._C()
        config = load_config(r)
        eq_(self.tempdir,config['NGSDATA'])

    def test_makes_config_to_specified_path(self, mock_yaml):
        from ngs_mapper.config import load_config
        mock_yaml.load.return_value = self.config
        os.mkdir('configdir')
        savepath = join('configdir', 'config.yaml')
        r = self._C(savepath)
        config = load_config(r)
        eq_(self.tempdir,config['NGSDATA'])
        eq_(r, savepath)

@patch('ngs_mapper.config.yaml')
@patch('__builtin__.open', MagicMock())
class TestGetConfigArgparse(Base):
    functionname = 'get_config_argparse'

    def test_outputs_version(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        with patch('ngs_mapper.config.sys') as msys:
            with patch('ngs_mapper.config.ngs_mapper') as mngs_mapper:
                mngs_mapper.__version__ = '1.1.1'
                msys.stdout = StringIO()
                r = self._C(['--version'])
                eq_('Version 1.1.1\n', msys.stdout.getvalue())

    def test_does_not_parse_help(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        r = self._C(['--help'])
        parser, args, config, configfile = r
        ok_('--help' in args, 'Parsed --help when it should not have')
    
    def test_returns_valid_argparse(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        from argparse import ArgumentParser
        r = self._C(['--foo', 'foo'])
        parser, args, config, configfile = r
        parser = ArgumentParser(parents=[parser])
        parser.add_argument('--foo')
        args = parser.parse_args(args)
        eq_(args.foo, 'foo')

    def test_returns_default_config(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        r = self._C([])
        parser, args, config, configfile = r
        eq_(self.tempdir, config['NGSDATA'])
        eq_(configfile, None)

    def test_returns_specified_config(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        with patch('ngs_mapper.config.load_config') as mock_load_config:
            r = self._C(['-c','/path/to/file.yaml'])
            parser, args, config, configfile = r
            mock_load_config.assert_called_once_with('/path/to/file.yaml')
            eq_('/path/to/file.yaml', configfile)

@patch('__builtin__.open')
@patch('ngs_mapper.config.yaml')
@patch('argparse.ArgumentParser.parse_args')
class TestMain(Base):
    functionname = 'main'

    def test_creates_config(self,mock_parse_args, mock_yaml, mock_open):
        curconfigpath = join(self.tempdir, 'my.yaml')
        mock_yaml.load.return_value = self.config
        parse_args = Mock()
        parse_args.save_to = curconfigpath
        mock_parse_args.return_value = parse_args

        import StringIO
        stdout = StringIO.StringIO()
        with patch('sys.stdout', stdout):
            self._C()
        # Make sure script outputs path to config created
        eq_(curconfigpath, stdout.getvalue().rstrip())
