from imports import *
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class Base( common.BaseClass ):
    modulepath = 'miseqpipeline.fqstats'

class TestReadAvgQual( Base ):
    functionname = 'read_avg_qual'

    def test_correct_avgq( self ):
        rrec = common.make_seqrec( 'ATGCATGC', [0,0,0,0,40,40,40,40] )
        r = self._C( rrec )
        eq_( 20.0, r )

    def test_rounds( self ):
        rrec = common.make_seqrec( 'ATGC', [0,0,5,5] )
        r = self._C( rrec )
        eq_( 3, r )

    def test_no_qual_values( self ):
        rrec = common.make_seqrec( '', [] )
        r = self._C( rrec )
        eq_( 0.0, r )

class TestBinValues( Base ):
    functionname = 'bin_values'

    def test_bins_correctly( self ):
        import random
        # List with unique random values
        values = list( set( [random.randint(1,100) for x in range(100)] ) )
        r, m = self._C( values * 10 )
        # Should be 10 for every value
        maxv = 0
        for v in values:
            eq_( 10, r[v] )
        eq_( 10, m, 'Maximum returned was not correct' )

class TestFqstats( Base ):
    functionname = 'fqstats'

    def test_expected_return( self ):
        acceptable_time = 0.5
        import time
        from miseqpipeline.fqstats import bin_values, read_avg_qual
        randseqs, mlen, mqual = common.random_seqs( 10000 )
        s = time.time()
        rlb, aqb, mrl, maq, mrc, mqc = self._C( randseqs )
        d = time.time() - s

        eq_( mlen, mrl, 'Did not compute max read length correctly {} != {}'.format(mlen,mrl) )
        eq_( mqual, maq, 'Did not compute max AvgQual correctly {} != {}'.format(mqual,maq) )
        rlbin, aa = bin_values( [len(r.seq._data) for r in randseqs] )
        aqbin, bb = bin_values( [read_avg_qual(r) for r in randseqs] )

        for k,v in rlbin.iteritems():
            eq_( v, rlb[k] )
        for k,v in aqbin.iteritems():
            eq_( v, aqb[k] )
        ok_( d < acceptable_time, 'Took {} seconds'.format(d) )
        eq_( aa, mrc, 'Did not set maximum read count correctly {} != {}'.format(100,mrc) )
        eq_( bb, mqc, 'Did not set maximum qual read count correctl {} != {}'.format(100,mqc) )

    def test_noseqs( self ):
        r = self._C( [] )
        eq_( ({},{},0,0,0,0), r )
