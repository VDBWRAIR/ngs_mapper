from os.path import *
import os
import sys

from mock import Mock, MagicMock, patch
from nose.tools import eq_, ok_, raises
from nose.plugins.attrib import attr

import common
import fixtures
from fixtures import THIS
