from os.path import *
import os
import sys
from cStringIO import StringIO
import tempfile
import shutil

from mock import Mock, MagicMock, patch
from nose.tools import eq_, ok_, raises
from nose.plugins.attrib import attr

import common
import fixtures
from fixtures import THIS
