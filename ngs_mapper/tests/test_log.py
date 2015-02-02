from nose.tools import eq_, ok_, raises
from mock import patch, Mock, MagicMock

import os
import sys
from os.path import *
import logging

import common

class LogBase(common.BaseClass):
    modulepath = 'ngs_mapper.log'

    def setUp( self ):
        super( LogBase, self ).setUp()
        self.filename = join('Projects','project','project.log')
        # Create intermediate directories up to path of project.log
        os.makedirs( dirname( self.filename ) )
        self.format = '%(asctime)-15s %(message)s'

class TestUnitGetConfig(LogBase):
    functionname = 'get_config'

    def is_valid_config( self, config ):
        import logging.config
        # Should throw exceptions if parsing fails
        logging.config.dictConfig( config )
    ivc = is_valid_config

    def test_valid_config( self ):
        r = self._C( self.filename, self.format )
        eq_( self.filename, r['handlers']['file']['filename'], "Did not set the filename for file handler" )
        eq_( self.format, r['formatters']['standard']['format'], "Did not set the format in standard formatter" )
        self.ivc( r )

class MockLoggingHandler(logging.Handler):
    def __init__( self, *args, **kwargs ):
        self.reset()
        logging.Handler.__init__( self, *args, **kwargs )

    def emit( self, record ):
        self.messages[record.levelname.lower()].append(record.getMessage())

    def reset( self ):
        self.messages = {
            'debug': [],
            'info': [],
            'warning': [],
            'error': [],
            'critical': []
        }

class TestUnitSetupLogger(LogBase):
    functionname = 'setup_logger'

    ''' Ensure logging happens as expected '''
    def setUp( self ):
        super( TestUnitSetupLogger, self ).setUp()
        from ngs_mapper.log import get_config
        self.config = get_config( self.filename, self.format )

    from StringIO import StringIO
    @patch('sys.stdout',new_callable=StringIO)
    def test_setup_logger( self, stdout ):
        r = self._C( __name__, self.config )
        r.info( "Test" )
        ok_( stdout.getvalue().endswith( "Test\n" ) )
        with open(self.filename) as fh:
            logs = fh.read()
            ok_( logs.endswith( "Test\n" ) )


import mock
import unittest2 as unittest

from os.path import *

from .. import log

class TestSetupLogger(unittest.TestCase):
    def setUp(self):
        self.filename = join('Projects','project','project.log')
        if not exists(dirname(self.filename)):
            os.makedirs(dirname(self.filename))
        self.format = '%(asctime)-15s %(message)s'
        self.config = log.get_config(self.filename, self.format)

    def test_python26_missing_dictconfig(self):
        with mock.patch.object(log.logging, 'config') as mlogconfig:
            mlogconfig.dictConfig.side_effect = ImportError
            log.setup_logger('foo', self.config)
