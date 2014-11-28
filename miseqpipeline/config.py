from os.path import exists, isdir, join
import os
import yaml
import pkg_resources
import argparse

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

def load_config(config_file):
    '''
    Loads a yaml config file from the given config_file path
     or from a stream
    Returns a dictionary
    '''
    if hasattr(config_file, 'read'):
        config_stream = config_file
    else:
        config_stream = open(config_file)
    config = yaml.load(config_file)
    verify_config(config)
    return config

def load_default_config():
    '''
    Loads the default config from pkg_resources
    Returns the config dictionary
    '''
    config_stream = pkg_resources.resource_stream(__name__, 'config.yaml')
    return load_config(config_stream)

def make_example_config(savepath=os.getcwd()):
    '''
    Load default config and dump to savepath
    '''
    if not exists(savepath):
        raise ValueError('{0} is not a valid path to save to'.format(savepath))
    if isdir(savepath):
        savepath = join(savepath, 'config.yaml')
    config = load_default_config()
    with open(savepath,'w') as fh:
        fh.write(
            yaml.dump(config, indent=4, default_flow_style=False)
        )
    return savepath

def get_config_argparse(argv):
    '''
    Setup an argparse instance for allowing user to specify the config file
    Then partially parse the command line args and return either the
    default config or the config the user specified

    Returns (config_parser, [rest of command line args], config)
    '''
    conf_parser = argparse.ArgumentParser(
        add_help=False
    )
    conf_parser.add_argument(
        '--config',
        '-c',
        dest='config',
        default=None,
        help='Path to config.yaml file'
    )

    args, rest = conf_parser.parse_known_args(argv)

    if args.config:
        config = load_config(args.config)
    else:
        config = load_default_config()

    return conf_parser, rest, config

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        description='Generate a config file that you can customize'
    )

    parser.add_argument(
        '--save-to',
        default=os.getcwd(),
        help='Where to generate the default config[Default: %s(default)]'
    )

    args = parser.parse_args()
