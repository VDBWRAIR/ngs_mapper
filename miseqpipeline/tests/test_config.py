from imports import *
from mock import mock_open

class Base(common.BaseClass):
    modulepath = 'miseqpipeline.config'

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

@patch('miseqpipeline.config.yaml')
class TestLoadConfigFile(Base):
    functionname = 'load_config'

    def test_loads_config_stream(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        config_stream = mock_open(read_data=self._create_yaml_from_config(self.config))()
        r = self._C(config_stream)
        eq_(self.tempdir, r['NGSDATA'])

    @patch('__builtin__.open', Mock())
    def test_loads_config_filepath(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        r = self._C('/path/to/file.yaml')
        eq_(self.tempdir, r['NGSDATA'])

    @patch('__builtin__.open', Mock())
    def test_invalid_config_raises_exception(self, mock_yaml):
        self.config['NGSDATA'] = '/missing/path'
        mock_yaml.load.return_value = self.config
        from miseqpipeline.config import InvalidConfigError
        assert_raises(InvalidConfigError, self._C, '/path/to/file.yaml')

class TestVerifyConfig(Base):
    functionname = 'verify_config'

    def test_bias_value_lt_1_raises_error(self):
        from miseqpipeline.config import InvalidConfigError
        self.config['base_caller']['bias']['default'] = 0
        assert_raises(InvalidConfigError, self._C, self.config)

    def test_ngsdata_not_set_raises_error(self):
        from miseqpipeline.config import InvalidConfigError
        self.config['NGSDATA'] = '/non/existant/path'
        assert_raises(InvalidConfigError, self._C, self.config)

    def test_valid_config_no_error(self):
        self.config['NGSDATA'] = self.tempdir
        self._C(self.config)

@patch('__builtin__.open', Mock())
@patch('pkg_resources.resource_stream')
@patch('miseqpipeline.config.yaml')
class TestLoadDefaultConfig(Base):
    functionname = 'load_default_config'

    def test_loads_from_pkg_resources(self, mock_yaml, mock_stream):
        mock_stream.return_value = mock_open(read_data=self._create_yaml_from_config(self.config))()
        mock_yaml.load.return_value = self.config
        r = self._C()
        eq_(self.tempdir, r['NGSDATA'])

@patch('miseqpipeline.config.yaml')
@patch('__builtin__.open', MagicMock())
@patch('pkg_resources.resource_string', Mock(return_value='/path/to/file.yaml'))
class TestMakeExampleConfig(Base):
    functionname = 'make_example_config'

    def test_invalid_savepath_raises_exception(self, mock_yaml):
        assert_raises(ValueError, self._C, '/non/existant/path.yaml')

    def test_makes_config_cwd_is_default(self, mock_yaml):
        from miseqpipeline.config import load_config
        mock_yaml.load.return_value = self.config
        r = self._C()
        config = load_config(r)
        eq_(self.tempdir,config['NGSDATA'])

    def test_makes_config_to_specified_path(self, mock_yaml):
        from miseqpipeline.config import load_config
        mock_yaml.load.return_value = self.config
        os.mkdir('configdir')
        savepath = join('configdir', 'config.yaml')
        r = self._C('configdir')
        config = load_config(r)
        eq_(self.tempdir,config['NGSDATA'])
        eq_(r, savepath)
