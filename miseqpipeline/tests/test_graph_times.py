from imports import *

class Base(BaseTester):
    modulepath = 'miseqpipeline.graph_times'

class TestHasNoTests(Base):
    def test_make_some_tests( self ):
        ok_( False, "No tests for this module yet" )
