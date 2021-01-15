from os.path import *
import os
import sys
from cStringIO import StringIO
import tempfile
import shutil
from glob import glob
import subprocess
import shlex
import re

from mock import Mock, MagicMock, patch, mock_open, call
from nose.tools import *
from nose.plugins.attrib import attr

from Bio import SeqIO

import common
import fixtures
from fixtures import THIS
from common import *
from . import tdir
