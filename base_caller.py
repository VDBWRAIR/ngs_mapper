#!/usr/bin/env python

#from stats_at_refpos import stats
from samtools import MPileupColumn, mpileup

from alphabet import iupac_amb
import sys
import argparse
import re
from Bio import SeqIO
from StringIO import StringIO
import vcf
from os.path import basename

# The header for the vcf
VCF_HEAD = '''##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=RC,Number=1,Type=Integer,Description="Reference Count">
##INFO=<ID=RAQ,Number=1,Type=Integer,Description="Reference Average Quality">
##INFO=<ID=PRC,Number=1,Type=Integer,Description="Percentage Reference Count">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternative Alternate Count">
##INFO=<ID=AAQ,Number=A,Type=Integer,Description="Alternate Average Quality">
##INFO=<ID=PAC,Number=A,Type=Integer,Description="Percentage Alternate Count">
##INFO=<ID=CBD,Number=1,Type=Integer,Description="Called Base Depth(Depth of only bases supporting CB">
##INFO=<ID=CB,Number=1,Type=Character,Description="Called Base">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}'''

def timeit( func ):
    def wrapper( *args, **kwargs ):
        import time; st = time.time()
        result = func( *args, **kwargs )
        print "{} took {} seconds to run".format(func, time.time() - st)
        return result
    return wrapper

def main( args ):
    generate_vcf(
        args.bamfile,
        args.reffile,
        args.regionstr,
        args.vcf_output_file,
        args.minbq,
        args.maxd,
        args.mind,
        args.minth,
        args.biasth,
        args.bias,
        VCF_HEAD.format(basename(args.bamfile)),
    )

def parse_args( args=sys.argv[1:] ):
    parser = argparse.ArgumentParser(
        description = 'Generates a VCF that has called bases in it which follow ' \
            'the WRAIR VDB SOP for calling bases',
        epilog = '''WRAIR VDB Base Calling SOP:
            Depth < 10:
                All bases with base quality < 25 get set to N
                Then the base is called on the percentage(See below)
            Depth > 10:
                All bases with base quality < 25 get removed
                Then the base is called on the percentage(See below)

            Calling on Percentage:
                Any base with >= 80% majority is called
                    - or -
                N is called if the depth was < 10 and N is > 20%
                    - or -
                The specific IUPAC ambiguious base is called for all bases over
                 20%
                    - or -
                The majority base is called
        '''
    )

    parser.add_argument(
        dest='bamfile',
        help='The bam file to generate the vcf for'
    )

    parser.add_argument(
        dest='reffile',
        help='The reference file(has to have an fai index already built)'
    )

    parser.add_argument(
        '-o',
        dest='vcf_output_file',
        default=None,
        help='Where to save the vcf'
    )

    parser.add_argument(
        '-r',
        '--regionstr',
        dest='regionstr',
        default=None,
        help='See the samtools documentation for how to specify a region string.' \
            'Essentially: \'refname\':START-STOP'
    )

    parser.add_argument(
        '-minbq',
        dest='minbq',
        default=25,
        type=int,
        help='The minimum base quality to be considered high quality[Default: 25]'
    )

    parser.add_argument(
        '-maxd',
        dest='maxd',
        default=100000,
        type=int,
        help='Maximum depth to use for the pileup[Default: 100000]'
    )

    parser.add_argument(
        '-mind',
        dest='mind',
        default=10,
        type=int,
        help='Minimum depth for base trimming. Below this depth low quality bases' \
            ' will be called N.[Default: 10]'
    )

    parser.add_argument(
        '-minth',
        dest='minth',
        default=0.8,
        type=float,
        help='Minimum fraction of all remaining bases after trimming/N calling that ' \
            'will trigger a base to be called.[Default: 0.8]'
    )

    parser.add_argument(
        '-biasth',
        dest='biasth',
        default=50,
        type=float,
        help='Minimum base quality threshold to bias towards. Will increase the amount of bases that have >= ' \
            'this value by a factor of what bias is set to[Default: 50]'
    )

    default=10
    parser.add_argument(
        '-bias',
        dest='bias',
        default=default,
        type=int,
        help='What factor to bias high quality bases by. Must be an integer >= 1[Default: {}]'.format(default)
    )

    args = parser.parse_args( args )
    if args.vcf_output_file is None:
        args.vcf_output_file = args.bamfile + '.vcf'

    return args

def mark_lq( stats, minbq, mind, refbase ):
    '''
        Labels all qualities < minbq as ?
        Goes through all keys in the stats dictionary that are not in ('depth','mqualsum','bqualsum')
        which should be keys that represent nucleotide bases. Those keys then point to a dictionary
        that contain 'mapq': [] and 'baseq': []

        Creating a new base called N or ? depending on the overall depth(N for < mind and ? for > mind)
        If the base is the reference base and < mind then the base will be preserved to bias the reference in low coverage areas

        @param stats - Stats dictionary returned from stats_at_refpos.stats
        @param minbq - The mininum base quality to determine if a quality should belong to N
        @param mind - The minimum depth threshold. If the depth is < this then lq will be labeled N otherwise
            they will be labeled ? for trimming purposes

        @returns stats dictionary with baseq subkey for each base in the dictionary.
    '''
    stats2 = {}
    stats2['depth'] = stats['depth']
    stats2['mqualsum'] = stats['mqualsum']
    stats2['bqualsum'] = stats['bqualsum']

    for base, quals in stats.iteritems():
        # Only interested in base stats in this loop
        if base not in ('depth','mqualsum','bqualsum'):
            # generates a list called bquals
            bquals = quals['baseq']
            # loop to examine the quality score and identifyes bases with a quality score less than the minbq of 25
            for q in bquals:
                k = base
                if q < minbq:
                    if stats2['depth'] < mind:
                        if base != refbase:
                            k = 'N'
                        else:
                            # Redundant but easier to see
                            k = base
                    else:
                        k = '?'
                # adds the N to the nucleotides (A C G T and N)
                if k not in stats2:
                    stats2[k] = {'baseq':[]}
                stats2[k]['baseq'].append( q )
    return stats2

def generate_vcf( bamfile, reffile, regionstr, vcf_output_file, minbq, maxd, mind=10, minth=0.8, biasth=50, bias=10, vcf_template=VCF_HEAD ):
    '''
        Generates a vcf file from a given vcf_template file

        @param bamfile - Path to bamfile
        @param reffile - Indexed reference path
        @param regionstr - samtools region string or None for all
        @param vcf_output_file - Where to write the output vcf
        @param minbq - Minumum base quality to determine if it should be turned into an N
        @param maxd - maximum depth to use
        @param vcf_template - VCF Header template(string)
        @param mind - Minimum depth to decide if low quality bases should be called N
        @param minth - Minimum percentage for a base to be called

        @returns path to vcf_output_file
    '''
    # All the references indexed by the seq.id(first string after the > in the file until the first space)
    refseqs = SeqIO.index( reffile, 'fasta' )
    # Our pretend file object that has vcf stuff in it
    vcf_head = StringIO( vcf_template )
    vcf_head.name = 'header.vcf'
    # Where to write the output file to
    if vcf_output_file is None:
        output_path = bamfile + '.vcf'
    else:
        output_path = vcf_output_file
    # The vcf writer object
    out_vcf = vcf.Writer( open( output_path, 'w' ), 
                            template = vcf.Reader( vcf_head )
    )

    # Get the iterator for an mpileupcal
    # Do not exclude any bases by setting minmq and minbq to 0 and maxdepth to 100000
    piles = mpileup( bamfile, regionstr, 0, 0, 100000 )
    # Loop through each pileup row
    # Last position stores the last position seen
    # Starts at 0 to show gap at beginning if needed
    lastpos = 0
    lastref = ''
    for pilestr in piles:
        # Generate the handy pileup column object
        col = MPileupColumn( pilestr )
        # Current position in alignment
        curpos = col.pos
        # The reference we are iterating on
        refseq = refseqs[col.ref].seq
        # Set lastref if we are on first iteration
        if lastref == '':
            lastref = col.ref
        # New reference so reset lastpos
        # and insert gaps at end
        elif lastref != col.ref:
            # Insert from lastposition of lastref to end of reference
            # len+1 so that it inserts single bases at the end
            for rec in blank_vcf_rows( lastref, refseq, len(refseq)+1, lastpos, '-' ):
                out_vcf.write_record( rec )
            lastref = col.ref
            lastpos = 0
        # Insert depth 0 regions
        for rec in blank_vcf_rows( col.ref, refseq, col.pos, lastpos, '-' ):
            out_vcf.write_record( rec )
        # Generate the vcf frow for that column
        row = generate_vcf_row( col, refseq, minbq, maxd, mind, minth, biasth, bias )
        # Write the record to the vcf file
        out_vcf.write_record( row )
        # Set last position seen
        lastpos = row.POS

    if lastref != '':
        # Insert from lastposition of lastref to end of reference
        # len+1 so that it inserts single bases at the end
        for rec in blank_vcf_rows( lastref, refseq, len(refseq)+1, lastpos, '-' ):
            out_vcf.write_record( rec )

    # Close the file
    out_vcf.close()

    return output_path

def blank_vcf_rows( refname, refseq, curpos, lastpos, call='-' ):
    '''
        Returns a list of blank vcf rows for all positions that are missing
        between curpos and lastpos.
        The blank rows should represent a gap in the alignment and be called whatever
        call is set too

        @param refname - Reference name to set _Record.CHROM
        @param refseq - Reference sequence to get reference base from
        @param curpos - Current position in the alignment(1 based)
        @param lastpos - Last position seen in the alignment(1 based)
        @param call - What to set the CB info field to(Default to -)

        @returns a list of vcf.model._Record objects filled out with the DP,RC,RAQ,PRC,CBD=0 and CB=call
    '''
    # Blank record list
    records = []
    # Only do records between curpos and lastpos
    for i in range( lastpos + 1, curpos ):
        # Reference base at current position
        rec = blank_vcf_row( refname, refseq, i, call )
        records.append( rec )

    return records

def blank_vcf_row( refname, refseq, pos, call='-' ):
    '''
        Generates a blank VCF row to insert for depths of 0

        @param refseq - Bio.seq.seq object representing the reference sequence
        @param pos - Reference position to get the reference base from
        @param call - What to set the CB info field to

        @returns a vcf.model._Record
    '''
    info = dict(
        DP=0,
        RC=0,
        RAQ=0,
        PRC=0,
        CB=call,
        CBD=0
    )
    record = vcf.model._Record( refname, pos, None, refseq[pos-1], '.', None, None, info, None, None )
    return record

def bias_hq( stats, biasth=50, bias=10 ):
    '''
        Biases high quality reads in stats so that they are more likely to be selected later on.
        Essentially duplicates quality scores(baseq values) by a factor of bias.
         Given the example:
            stats = { 'A': {'baseq':[40,40,40]}, 'C': {'baseq':[50,50,50,50,50,50,50]} }
         Without any biasing the called base would likely be called an M as there is no majority of C's
          since A is 30% and C is 70%(Assuming 0.8 minth which requires 80% to be called)
         With bias_hq( stats, 50, 2 ) you would end up with 3 A's and 14 C's or 18% A and 82% C so C
            would be more correctly called for the consensus

        bias will round the length of the created list to the nearest integer so if the length of a
        baseq is say 10 and the bias is 1.25 the resulting baseq list will be 13 in length(addition of 3 baseq)

        @param stats - Should most likely be a stats2 dictionary although stats would work as well
        @param biasth - What quality value(>=) should be considered to be bias towards
        @param bias - How much to bias aka, how much to multiply the # of quals >= biasth(has to be int >= 1)

        @returns stats2 formatted dictionary with all baseq lists appended to with the bias amount
    '''
    if bias < 1 or int(bias) != bias:
        raise ValueError( "bias was set to {} which is less than 1. Cannot bias on a factor < 1".format(bias) )

    # Because we are adding onto the existing list this needs to be one less
    bias = int(bias) - 1

    stats2 = {'depth': 0}

    for k, v in stats.iteritems():
        # Skip non base items
        if k in ('depth','mqualsum','bqualsum'):
            if k != 'depth':
                stats2[k] = v
            continue
        stats2[k] = {}
        bquals = v['baseq']
        # New list of all baseq >= biasth biased by a factor of bias
        duplist = [d for d in bquals if d >= biasth] * bias
        stats2[k]['baseq'] = v['baseq'] + duplist
        stats2[k]['mapq'] = v.get( 'mapq', [] )
        stats2['depth'] += len( stats2[k]['baseq'] )
    return stats2

def pile_stats( mpileupcol, refbase, minbq, mind, biasth, bias ):
    '''
        Returns the modified statistics from an mpileupcol that are suitable for base calling

        All parameters are identical to generate_vcf_row

        @returns a stats2 dictionary that is modified by biasing reference bases and high quality bases
    '''
    # This call needs to be replaced after the mpileupcol parameter is put into place
    # then stats_for_col( mpileup
    s = mpileupcol.base_stats()
    # Bias high quality first as it may change the behavior of mark_lq as the depth may
    # increase above the mind threshold
    stats2 = bias_hq( s, biasth, bias )
    stats2 = mark_lq( stats2, minbq, mind, refbase )

    # Update stats2 so that it does not include low quality bases since we
    # are equal to or above the min depth
    if stats2['depth'] >= mind:
        if '?' in stats2:
            stats2['depth'] -= len(stats2['?']['baseq'])
            del stats2['?']

    return stats2

def generate_vcf_row( mpileupcol, refseq, minbq, maxd, mind=10, minth=0.8, biasth=50, bias=10 ):
    '''
        Generates a vcf row and returns it as a string

        @param mpileupcol - samtools.MpileupColumn
        @param refseq - Bio.Seq.Seq object representing the reference sequence
        @param minbq - minimum base quality to be considered or turned into an N
        @param maxd - Maximum depth for pileup
        @param mind - Minimum depth decides if low quality bases are N's or if they are removed
        @param minth - Minimum percentage to call a base(unless no bases have > minth then the maximum pct base would be called
        @param biasth - Any base with quality above this will be biased by a bias factor
        @param bias - For every base >= biasth add bias more of those bases

        @returns a vcf.model._Record
    '''
    from collections import OrderedDict
    # The base position should be the same as the second item in the parsed region string
    start = mpileupcol.pos
    # Get the reference base from the reference sequence
    # Python is 0-index, biology is 1 index
    rb = refseq[start-1].upper()

    stats2 = pile_stats( mpileupcol, rb, minbq, mind, biasth, bias )

    # info needs to contrain the depth, ref count #, % ref count, Ave ref qual, alt count #, % ref count, Ave alf qual
    # Holds the info dictionary in order
    info = {}
    # Alternate info
    alt_info = info_stats( stats2, rb )
    alt_bases = alt_info['bases']
    # Only merge if there is something to merge
    if alt_info['bases'] != []:
        del alt_info['bases']
        # Merge existing info, with alt_info
        info['AC'] = alt_info['AC']
        info['AAQ'] = alt_info['AAQ']
        info['PAC'] = alt_info['PAC']

    # Total Depth at this position
    info['DP'] = stats2['depth']
    # Quick alias to reference statistics
    # Maybe reference base isn't in stats
    if rb in stats2:
        refstats = stats2[rb]
        # Reference Count is length of the base qualities list
        info['RC'] = len(refstats['baseq'])
        # Reference Average Quality is the sum of base qualities / length
        info['RAQ'] = int( round( sum(refstats['baseq']) / float(len(refstats['baseq'])), 0 ) )
        # Percentage Reference Count is len of qualities / depth
        info['PRC'] = int( round( (100.0 * len(refstats['baseq'])) / float(stats2['depth']), 0 ) )
    else:
        info['RC'] = 0
        info['RAQ'] = 0
        info['PRC'] = 0

    # Get Called Base Info
    cb, cbd = caller( stats2, minbq, maxd, mind, minth )

    # Set the Called base
    info['CB'] = cb
    info['CBD'] = cbd

    if not alt_bases:
        alt_bases = '.'

    # need to record each line of the vcf file.
    record = vcf.model._Record( mpileupcol.ref, start, None, rb, alt_bases, None, None, info, None, None )

    return record

def caller( stats2, minbq, maxd, mind=10, minth=0.8 ):
    '''
        Calls a given base at refstr inside of bamfile. At this time refstr has to be a single
        base position('refname':N-N). The base is determined by first labeling all bases less than minbq as N and then
        determining if the depth is < mind or >= mind.
        If < and the % of N is > minth then call it an N as it is the majority.
        If >= 10 then remove all N

        The final stage is to call call_on_pct with the remaining statistics

        @param stats - stats dictionary generated from stats_at_refpos.stats
        @param minbq - Minimum base quality to determine if it is low quality or not to call it an N
        @param maxd - Maximum depth in the pileup to allow
        @param mind - Minimum depth threshold
        @param minth - Minimum percentage for an base to be called non ambigious

        @returns one of the items from the set( 'ATGCMRWSYKVHDBN' )
    '''
    # Else we will determine if the N's are the majority now
    # defines if the base is an N
    if '?' in stats2:
        nlen = len(stats2['?']['baseq'])
        np = nlen/(stats2['depth']*1.0)
        if np > (1-minth):
            return ('N', nlen)
    return call_on_pct(stats2, minth)

def call_on_pct( stats2, minth=0.8 ):
    '''
        Calls a base from the given stats dictionary if it is the majority. A majority base is
        any base where it is in %total >= minth.

        Base quality is not used to make the determination, but instead the number of quality scores
        in the baseq list is used as it depicts the depth of that base. The quality scores are not used
        as the stats dictionary should already be run through the label_n function.

        @param stats2 - Stats dictionary returned from mark_lq or stats_at_refpos.stats
        @param minth - minimum percentage that a base needs to be present in order to be called non-ambiguous

        @returns the called base based on the percentages in the given stats and the depth for the called base
    '''
    # Empty stats indicates zero depth
    if not stats2:
        return ('-',0)
    nt_list = ''
    count = 0
    for base, quals in stats2.iteritems():
        # Only interested in base stats in this loop
        if base not in ('depth','mqualsum','bqualsum'):
            # Quick alias for quals['baseq']
            bquals = quals['baseq']    
            # Percentage of current base compared to total depth
            np_2 = len(bquals)/(stats2['depth']*1.0)
            # fix for proper calculations with float
            # If basepercent is greater than minimum threashold
            if np_2 > round((1-minth),1):
                nt_list += base
                count += len(bquals)
    dnalist = sorted(nt_list)
    try:
        return (iupac_amb(dnalist), count)
    except ValueError as e:
        # Only depth from bases for the N
        if dnalist:
            ndepth = sum([len(stats2[b]['baseq']) for b in dnalist])
        else:
            ndepth = stats2.get( 'depth', 0 )
        return ('N', ndepth)

def info_stats( stats, rb):
    '''
        Returns a dictionary that can be used to fill in the info field of a vcf record
        The returned dictionary will contain the AC(Alternate Count), AAQ(Alternate Average Quality)
            PAC(Percentage Average Count) and bases.

        AC should be a list of integers
        AAQ should be a list of integers
        PAC should be a list of integers rounded to the nearest integer
        bases should be the ordered list of bases for each item in AC,AAQ, PAC such that
        zip( bases, AC/AAQ/PAC ) would work as expected
    
        @param stats - stats dictionary(should accept either stats or stats2)
        @param rb - Reference base to ignore in the outputted dictionary

        @returns info dictionary with AC,AAQ, PAC and bases keys filled out
    '''
    info = {
        'AC' : [], 
        'AAQ' : [],
        'PAC' : [],
        'bases' : []
    }

     # the alternitive base = ab
    alt_nt = []
    for base, quals in stats.iteritems():
        if base  not in ('depth','mqualsum','bqualsum',rb):
            # identify the alternitive bases in stats 2        
            # data for the alternitive count
            info['AC'].append(len(quals['baseq']))
            # data for the alternitive avarage quality
            info['AAQ'].append(int(round((sum(quals['baseq'])*1.0)/len(quals['baseq']))))
            # data for the percentage reference count
            info['PAC'].append(int(round((len(quals['baseq'])*100.0)/(stats['depth']))))
            # base data
            info['bases'].append(base)
                                     
    return info
 
if __name__ == '__main__':
    main( parse_args() )
