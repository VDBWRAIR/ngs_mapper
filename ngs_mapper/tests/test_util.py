import os.path
from os.path import *

import unittest
import mock

from .. import util

@mock.patch.object(util, 'os')
class TestBuildDatafiles(unittest.TestCase):
    def setUp(self):
        # For os.walk
        self.walk = [
            ['', ['dir1','dir2'], ['1.txt','2.txt']],
            ['dir1', [], ['3.txt','4.txt']],
            ['dir2', ['dir3'], ['5.txt', '6.txt']],
            ['dir2/dir3', [], ['7.txt']],
        ]

    def _oswalk(self, path):
        for rdf in self.walk:
            rdf[0] = join(path, rdf[0])
            yield tuple(rdf)

    def test_returns_correct_manifest(self, mos):
        mos.walk.side_effect = self._oswalk
        mos.path = os.path
        r = util.build_datafiles('prefix', 'root')
        self.assertEqual(
            [
                ('prefix', ['root/1.txt','root/2.txt']),
                ('prefix/dir1', ['root/dir1/3.txt','root/dir1/4.txt']),
                ('prefix/dir2', ['root/dir2/5.txt','root/dir2/6.txt']),
                ('prefix/dir2/dir3', ['root/dir2/dir3/7.txt']),
            ],
            r
        )
