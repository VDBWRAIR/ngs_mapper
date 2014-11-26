from os.path import exists
import os
import yaml
import pkg_resources

# Raised when invalid config is loaded
class InvalidConfigError(Exception): pass

def verify_config(config):
    '''
    Verifies a config dictionary to ensure that specific values
    exist and are valid

    Raises InvalidConfigError if any errors are found
    '''
    if not exists(config['NGSDATA']):
        raise InvalidConfigError(
            "{0} is not a valid NGSDATA path".format(config['NGSDATA'])
        )

    if config['base_caller']['bias']['default'] < 1:
        raise InvalidConfigError(
            "base_caller:bias needs to be an integer >= 1".format(
                config['base_caller']['bias']['default']
            )
        )

def load_config_file(config_file):
    '''
    Loads a yaml config file from the given config_file path
    Returns a dictionary
    '''
    config = yaml.load(open(config_file))
    verify_config(config)
    return config

def load_default_config():
    '''
    Loads the default config from pkg_resources
    Returns the config dictionary
    '''
    config_path = pkg_resources.resource_string(__name__, 'config.yaml')
    return load_config_file(config_path)

def make_example_config(savepath=os.getcwd()):
    '''
    Load default config and dump to savepath
    '''
    config = load_default_config()
    with open(savepath,'w') as fh:
        fh.write(yaml.dump(config))
    return savepath
