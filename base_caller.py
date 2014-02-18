#!/usr/bin/env python

from stats_at_refpos import stats
from alphabet import iupac_amb
import sys
import argparse
import re
import pysam
from Bio import SeqIO
from StringIO import StringIO

# The header for the vcf
VCF_HEAD = '''##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=RC,Number=1,Type=Integer,Description="Reference Count">
##INFO=<ID=ARQ,Number=1,Type=Float,Description="Average Reference Quality">
##INFO=<ID=PRC,Number=1,Type=Integer,Description="Percentage Reference Count">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternative Alternate Count">
##INFO=<ID=AAQ,Number=A,Type=Float,Description="Alternate Average Quality">
##INFO=<ID=PAC,Number=A,Type=Integer,Description="Percentage Alternate Count">
##INFO=<ID=CBD,Number=1,Type=Integer,Description="Called Base Depth(Depth of only bases supporting CB">
##INFO=<ID=CB,Number=1,Type=Character,Description="Called Base">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'''

def main( args ):
    generate_vcf(
        args.bamfile,
        args.reffile,
        args.regionstr,
        args.vcf_output_file,
        args.minbq,
        args.maxd,
        VCF_HEAD,
        args.mind,
        args.minth
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

    args = parser.parse_args( args )
    if args.vcf_output_file is None:
        args.vcf_output_file = args.bamfile + '.vcf'

    return args

def label_N( stats, minbq ):
    '''
        Labels all qualities < minbq as N
        Goes through all keys in the stats dictionary that are not in ('depth','mqualsum','bqualsum')
        which should be keys that represent nucleotide bases. Those keys then point to a dictionary
        that contain 'mapq': [] and 'baseq': []

        This script loops through baseq and for any value less than minbq it will remove that quality score
        and place it in a new base key for 'N'.

        Essentially creating a new base N with all < minbq base qualities.

        This means that the mapq scores are all removed in the returned dictionary

        @param stats - Stats dictionary returned from stats_at_refpos.stats
        @param minbq - The mininum base quality to determine if a quality should belong to N

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
                    k = 'N'
                # adds the N to the nucleotides (A C G T and N)
                if k not in stats2:
                    stats2[k] = {'baseq':[]}
                stats2[k] ['baseq'].append( q )
    return stats2

def generate_vcf( bamfile, reffile, regionstr, vcf_output_file, minbq, maxd, vcf_template, mind=10, minth=0.8 ):
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
    import vcf
    refseqs = SeqIO.index( reffile, 'fasta' )
    # Do all references
    if regionstr is None:
        regions = [(r,1,len(seq.seq)) for r,seq in refseqs]
    else: # Do only single reference
        region = parse_regionstring( regionstr )
        regions = [region]
    # Our pretend file object that has vcf stuff in it
    vcf_head = StringIO( VCF_HEAD )
    vcf_head.name = 'header.vcf'
    # The vcf writer object
    out_vcf = vcf.Writer( open( vcf_output_file, 'w' ), 
                            template = vcf.Reader( vcf_head )
    )
    # Loop through all the positions in regionstr
    for ref, start, end in regions:
        for i in range( start, end+1 ):
            # Generate the new region string
            regionstr = region[0] + ':' + str(i) + '-' + str(i)
            refseq = refseqs[region[0]]
            # Probably inefficient that we are sending in the bamfile that has to be opened over and
            # over, but for now we will do that
            row = generate_vcf_row( bamfile, regionstr, refseq, minbq, maxd, mind, minth )
            out_vcf.write( row )

    return vcf_output_file

# Exception for when invalid region strings are given
class InvalidRegionString(Exception): pass

def parse_regionstring( regionstr ):
    '''
        Parses a region string into a 3 item tuple and checks it for errors

        @param regionstr - samtools region string format

        @returns (ref, start, stop)
        @raises InvalidRegionString
    '''
    m = re.match( '(\S+):(\d+)-(\d+)', regionstr )
    if not m:
        raise InvalidRegionString( "{} is not a valid regionstring".format(regionstr) )

    region = (m.group(1), int(m.group(2)), int(m.group(3)))
    if region[1] > region[2]:
        raise InvalidRegionString( "Start cannot be > stop in a region string: {}".format(regionstr) )

    return region

def generate_vcf_row( bam, regionstr, refseq, minbq, maxd, mind=10, minth=0.8 ):
    '''
        Generates a vcf row and returns it as a string

        @param bam - pysam.Samfile object
        @param regionstr - samtools regionstring to look at for a specific base(aka refname:N-N)
        @param refseq - Bio.Seq.Seq object representing the reference sequence
        @param minbq - minimum base quality to be considered
        @param maxd - Maximum depth for pileup
        @param mind - Minimum depth decides if low quality bases are N's or if they are removed
        @param minth - Minimum percentage to call a base(unless no bases have > minth then the maximum pct base would be called

        @returns a vcf.model._Record
    '''



   stats2 = label_N( stats (bamfile, regionstr, minmq, minbq, maxd ):
 



    # info needs to contrain the depth, ref count #, % ref count, Ave ref qual, alt count #, % ref count, Ave alf qual
    info = {
        'DP': 0,
        'RC': 0,
        'RAQ': 0,
        'PRC': 0,
        'AC': 0,
        'AAQ': 0,
        'PAC': 0,
        'CB': 0,
        'CBQ': 0,
    }


    # retrieve the info directly from the stats2 dictionary
    info['DP'] = stats2['depth']
    
    # find the base on the reference file - use this by the parse command ensure selection of the base    
    r = parse_regionstr(regionstr)
    rb = refseq[r[1]]


    # data for the depth
    info['DP'] = stats['depth']

    # data for the reference count
    info['RC'] = len(stats[rb])

    # data for the reference average quality
    info['RAQ'] = sum(stats[rb]['baseq'])/len(stats[rb]['baseq'])

    # data for the percentage reference count
    info['PRC'] = len(stats[rb])/(stats['depth']*1.0)



    # need to generate the stats2 information for the alternitive nucleotide - 
    # need to remove the data that is linked to the reference nucleotide




 

    


    # need to record each line of the vcf file.
    record = vcf._Record( 


        
    


    # call the nucleotides
    #s stats(bamfile, regionstr, 1, 1, maxd)
    pass

def caller( stats, minbq, maxd, mind=10, minth=0.8 ):
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
    # Label low quality as N
    stats2 = label_N(stats, minbq)
    # if the depth is >= min depth then just remove the N's
    if stats2['depth'] >= mind:
        if 'N' in stats2:        
            del stats2['N']
    else:
        # Else we will determine if the N's are the majority now
        # defines if the base is an N
        if 'N' in stats2:
            np = len(stats2['N']['baseq'])/(stats2['depth']*1.0)
            if np > (1-minth):
                return 'N'
    return call_on_pct(stats2, minth)

def call_on_pct( stats, minth=0.8 ):
    '''
        Calls a base from the given stats dictionary if it is the majority. A majority base is
        any base where it is in %total >= minth.

        Base quality is not used to make the determination, but instead the number of quality scores
        in the baseq list is used as it depicts the depth of that base. The quality scores are not used
        as the stats dictionary should already be run through the label_n function.

        @param stats2 - Stats dictionary returned from label_N or stats_at_refpos.stats
        @param minth - minimum percentage that a base needs to be present in order to be called non-ambiguous

        @returns the called base based on the percentages in the given stats
    '''
    nt_list = []
    for base, quals in stats.iteritems():
        # Only interested in base stats in this loop
        if base not in ('depth','mqualsum','bqualsum'):
            # generates a list called bquals
            bquals = quals['baseq']    
            np_2 = len(bquals)/(stats['depth']*1.0)
            if np_2 > round((1-minth),1): # fix for proper calculations with float
                nt_list.append( base )
    dnalist = ''.join(sorted(nt_list))
    return iupac_amb(dnalist)
    

def info_alt_stats( stats2, rb)


    info = {
        'AC' : [], 
        'AAQ' : [],
        'PAC' : [],
        'bases' : []
        }

     # the alternitive base = ab
    alt_nt = []
    for base, quals in stats.intritems():
        if not in ('depth','rb')
            ab = base  


    # identify the alternitive bases in stats 2        

            # data for the alternitive count
            info['AC'].append(len(quals['baseq']))


            # data for the alternitive avarage quality
            info['AAQ'].append(sum(quals['baseq'])/quals['baseq']))

        
            # data for the percentage reference count
            info['PAC'].append(len(quals)/(quals['depth']*1.0))
 

            # base data
            info['bases'].append(base)
                                     
    return info



 
if __name__ == '__main__':
    main( parse_args() )

