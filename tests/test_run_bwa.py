from nose.tools import eq_, raises
from nose.plugins.attrib import attr
from mock import MagicMock, patch

from common import BaseClass

import os

class Base(BaseClass):
    pass

class TestFunctional(Base):
    pass

class TestIntegration(Base):
    pass
