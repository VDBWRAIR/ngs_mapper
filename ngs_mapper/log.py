import logconfig
import logging
import logging.config

def setup_logger( name, config ):
    '''
        Returns a logging.getLogger instance with name and configured with dictConfig using config

        @param name - Name to give the logger(passed to logging.getLogger
        @param config - Valid logging.config.dictConfig object(easiest to use get_config to get this)

        @returns a logging.getLogger instance that is configured
    '''
    # Add name as a logger
    config['loggers'][name] = {
        'level': 'DEBUG',
        'handlers': config['handlers'].keys()
    }
    #logging.config.dictConfig( config )
    logconfig.from_dict(config)
    logging.getLogger("sh").setLevel(logging.WARNING)
    log = logging.getLogger( name )
    return log

def get_config( filename='pipeline.log', format='%(asctime)-15s -- %(levelname)s -- %(name)-15s %(message)s' ):
    '''
        Returns a standardized logging.config.dictConfig parsable dictionary
        with values set appropriate for this pipeline

        @param filename - The filename to log DEBUG and above to(probably inside of a project directory)
        @param format - Format value

        @returns a valid logging config dictionary
    '''
    log_config = {
        'version': 1,
        'handlers': {
            'console': {
                'class': 'logging.StreamHandler',
                'formatter': 'standard',
                'level': 'INFO',
                'stream': 'ext://sys.stdout'
            },
            'file': {
                'class': 'logging.FileHandler',
                'formatter': 'standard',
                'level': 'DEBUG',
                'filename': filename
            }
        },
        'formatters': {
            'standard': {
                'format': format
            }
        },
        'loggers': {
            'root': {
                'level': 'DEBUG',
                'handlers': ['console','file']
            }
        }
    }

    return log_config
