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
@patch('__builtin__.open', Mock())
class TestLoadConfigFile(Base):
    functionname = 'load_config_file'

    def test_loads_config(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        r = self._C('/path/to/file.yaml')
        eq_(self.tempdir, r['NGSDATA'])

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

@patch('miseqpipeline.config.yaml')
@patch('__builtin__.open', Mock())
@patch('pkg_resources.resource_string', Mock(return_value='/path/to/file.yaml'))
class TestLoadDefaultConfig(Base):
    functionname = 'load_default_config'

    def test_loads_from_pkg_resources(self, mock_yaml):
        mock_yaml.load.return_value = self.config
        r = self._C()
        eq_(self.tempdir, r['NGSDATA'])

@patch('miseqpipeline.config.yaml')
@patch('__builtin__.open', MagicMock())
@patch('pkg_resources.resource_string', Mock(return_value='/path/to/file.yaml'))
class TestMakeExampleConfig(Base):
    functionname = 'make_example_config'

    def test_makes_config_cwd_is_default(self, mock_yaml):
        from miseqpipeline.config import load_config_file
        mock_yaml.load.return_value = self.config
        r = self._C()
        config = load_config_file(r)
        eq_(self.tempdir,config['NGSDATA'])

    def test_makes_config_to_specified_path(self, mock_yaml):
        from miseqpipeline.config import load_config_file
        mock_yaml.load.return_value = self.config
        os.mkdir('configdir')
        r = self._C('configdir')
        config = load_config_file(r)
        eq_(self.tempdir,config['NGSDATA'])
