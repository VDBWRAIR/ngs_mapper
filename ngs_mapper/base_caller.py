from ngs_mapper.samtools import MPileupColumn, mpileup, parse_regionstring
from ngs_mapper.alphabet import iupac_amb

import sys
import argparse
import re
from StringIO import StringIO
from os.path import basename
import os
import multiprocessing
import time

import vcf
from Bio import SeqIO

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
##INFO=<ID=HPOLY,Number=0,Type=Flag,Description="Is a homopolymer">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}'''

def timeit(func):
    def wrapper(*args, **kwargs):
        import time; st = time.time()
        result = func(*args, **kwargs)
        print "{} took {} seconds to run".format(func, time.time() - st)
        return result
    return wrapper

def main():
    args = parse_args()
    if args.regionstr is not None:
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
            True
       )
    else:
        generate_vcf_multithreaded(
                args.bamfile,
                args.reffile,
                args.vcf_output_file,
                args.minbq,
                args.maxd,
                args.mind,
                args.minth,
                args.biasth,
                args.bias,
                args.threads,
                VCF_HEAD.format(basename(args.bamfile)),
       )

def generate_vcf_multithreaded(bamfile, reffile, vcf_output_file, minbq, maxd, mind, minth, biasth, bias, threads, vcfhead=VCF_HEAD):
    '''
    Generate vcf for each ref and split each ref into pieces
    '''
    # Generate name if not given
    if vcf_output_file is None:
        vcf_output_file = bamfile + '.vcf'
    # Store the temp file names to concat together later
    map_args = []
    # Temporary name suffix because tmpfile is too good of an idea
    i = 0

    procs = []
    for refrecord in SeqIO.parse(reffile,'fasta'):
        reflen = len(refrecord.seq)+1

        # Will break reference into chunks depending on how many threads
        chunksize = reflen/threads
        
        # Make regionstr to match chunks
        for start in range(1, reflen, chunksize):
            vcf_tmp_filename = "{0}.{1}".format(vcf_output_file,i)
            end = start + chunksize - 1
            regionstr = '{0}:{1}-{2}'.format(refrecord.id,start,end)
            
            complete_ref = False
            if end == reflen:
                complete_ref = True

            args = (bamfile, reffile, regionstr, vcf_tmp_filename, minbq, maxd, mind, minth, biasth, bias, vcfhead)#, complete_ref)
            map_args.append(args)
            # Create and start a new process
            p = multiprocessing.Process(target=generate_vcf, args=args)
            p.start()
            procs.append(p)
            #generate_vcf(*args)
            i += 1

        # Wait for all processes to finish
        for p in procs:
            p.join()

    with open(vcf_output_file, 'w') as fho:
        # Write the head
        fho.write(vcfhead + '\n')
        # Cat all tmpfiles and remove them
        for args in map_args:
            f = args[3]
            # Write tmp file contents to master file
            with open(f) as fhr:
                # Burn off header
                for line in VCF_HEAD.splitlines():
                    fhr.readline()
                fho.write(fhr.read())
                # Remove temp file
                os.unlink(f)
    return vcf_output_file

def parse_args(args=sys.argv[1:]):
    from ngs_mapper import config
    conf_parser, args, config, configfile = config.get_config_argparse(args)
    defaults = config['base_caller']

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
        ''',
        parents=[conf_parser]
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
        dest='vcf_output_file',
        help='Where to save the vcf'
   )

    parser.add_argument(
        '-r',
        '--regionstr',
        dest='regionstr',
        default=defaults['regionstr']['default'],
        help=defaults['regionstr']['help']
   )

    parser.add_argument(
        '-minbq',
        dest='minbq',
        default=defaults['minbq']['default'],
        type=int,
        help=defaults['minbq']['help']
   )

    parser.add_argument(
        '-maxd',
        dest='maxd',
        default=defaults['maxd']['default'],
        type=int,
        help=defaults['maxd']['help']
   )

    parser.add_argument(
        '-mind',
        dest='mind',
        default=defaults['mind']['default'],
        type=int,
        help=defaults['mind']['help']
   )

    parser.add_argument(
        '-minth',
        dest='minth',
        default=defaults['minth']['default'],
        type=float,
        help=defaults['minth']['help']
   )

    parser.add_argument(
        '-biasth',
        dest='biasth',
        default=defaults['biasth']['default'],
        type=float,
        help=defaults['biasth']['help']
   )

    parser.add_argument(
        '-bias',
        dest='bias',
        default=defaults['bias']['default'],
        type=int,
        help=defaults['bias']['help']
   )

    parser.add_argument(
        '--threads',
        default=defaults['threads']['default'],
        type=int,
        help=defaults['threads']['help']
   )

    args = parser.parse_args(args)
    if args.vcf_output_file is None:
        args.vcf_output_file = args.bamfile + '.vcf'

    return args

def mark_lq(stats, minbq, mind, refbase):
    '''
    Goes through all keys in the stats dictionary that are not in ('depth','mqualsum','bqualsum')
    which should be keys that represent nucleotide bases. Those keys then point to a dictionary
    that contain 'mapq': [] and 'baseq': []
    Creating a new base called N or ? depending on the overall depth

        - N for < mind & not refbase
        - Base for < mind & refbase
        - ? for > mind

    If the base is the reference base and < mind then the base will be preserved to bias the reference in low coverage areas

    :param str stats: Stats dictionary returned from stats_at_refpos.stats
    :param str minbq: The mininum base quality to determine if a quality should belong to N
    :param str mind: The minimum depth threshold. If the depth is < this then lq will be labeled N otherwise they will be labeled ? for trimming purposes

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
            # loop to examine the quality score and identifyes bases with a quality score less than the minbq
            for q in bquals:
                # Copy base so we can modify if needed without destroying original base
                k = base
                if q < minbq:
                    # Determine base to use since this base is low quality
                    if stats2['depth'] < mind:
                        if base != refbase:
                            # N since low qual and low depth
                            k = 'N'
                        else:
                            # Bias reference base
                            k = base
                    else:
                        # Base is unknown
                        k = '?'
                # adds the N to the nucleotides (A C G T and N)
                if k not in stats2:
                    stats2[k] = {'baseq':[]}
                stats2[k]['baseq'].append(q)
    return stats2

def hpoly_list(refseqs, minlength=3):
    '''
    Identify all homopolymer regions inside of each sequence in refseqs

    :param str refseqs: Bio.SeqIO.index'd fasta
    '''
    hpolys = {}
    p = r'([ATGC])\1{'+str(minlength-1)+',}'
    for seq in refseqs:
        matches = re.finditer(p, str(refseqs[seq].seq))
        hpolys[seq] = [(m.group(0),m.start()+1,m.end()) for m in matches]
    return hpolys

def is_hpoly(hpolylist, seqid, curpos):
    '''
    Identifies if a position is contained inside of a homopolymer
    '''
    l = hpolylist[seqid]
    for polys in l:
        if curpos >= polys[1] and curpos <= polys[2]:
            return True
    return False

def generate_vcf(bamfile, reffile, regionstr, vcf_output_file, minbq, maxd, mind=10, minth=0.8, biasth=50, bias=10, vcf_template=VCF_HEAD, complete_ref=False):
    '''
    Generates a vcf file from a given vcf_template file

    30000 bases ~= .05 seconds
    300000 bases ~= .5 seconds
    3000000 bases ~= 5 seconds
    30000000 bases ~= 50 seconds

    :param str bamfile: Path to bamfile
    :param str reffile: Indexed reference path
    :param str regionstr: samtools region string
    :param str vcf_output_file: Where to write the output vcf
    :param int minbq: Minumum base quality to determine if it should be turned into an N
    :param int maxd: maximum depth to use
    :param int mind: Minimum depth to decide if low quality bases should be called N
    :param float minth: Minimum percentage for a base to be called
    :param int biasth: What quality value(>=) should be considered to be bias towards
    :param str bias: For every base >= biasth add bias more of those bases
    :param str vcf_template: VCF Header template(string)
    :param bool complete_ref: If True, then complete all the way to the end position in regionstr

    @returns path to vcf_output_file
    '''
    #print regionstr
    # All the references indexed by the seq.id(first string after the > in the file until the first space)
    refseqs = SeqIO.index(reffile, 'fasta')
    # Homopolymers for references
    hpolys = hpoly_list(refseqs, 3)
    # Our pretend file object that has vcf stuff in it
    vcf_head = StringIO(vcf_template)
    vcf_head.name = 'header.vcf'
    # Where to write the output file to
    if vcf_output_file is None:
        output_path = bamfile + '.vcf'
    else:
        output_path = vcf_output_file
    # The vcf writer object
    fh = open(output_path, 'w')
    out_vcf = vcf.Writer(
        fh,
        template=vcf.Reader(vcf_head)
    )


    # Get the iterator for an mpileupcal
    # Do not exclude any bases by setting minmq and minbq to 0 and maxdepth to 100000
    piles = mpileup(bamfile, regionstr, 0, 0, 100000)

    # Parse the region string for later
    parsed_regionstr = parse_regionstring(regionstr)
    # Get the reference name to work with
    refname = parsed_regionstr[0]
    refseq = refseqs[refname].seq
    # The end of the ref may be restricted via regionstr
    # Lets user specify region start less than 1
    refstart = max(parsed_regionstr[1], 1)
    # Set refstart to 0 if regionstring start is 1
    # so that the first base gets inserted with blank_vcf_rows
    if refstart == 1:
        refstart = 0
    # Lets user specify region end past length of region
    refend = min(parsed_regionstr[2], len(refseq))
    # Last position stores the last position seen
    lastpos = refstart

    # Loop through each pileup row
    for pilestr in piles:
        # Generate the handy pileup column object
        col = MPileupColumn(pilestr)
        # Current position in alignment
        curpos = col.pos
        # The reference we are iterating on
        for rec in blank_vcf_rows(col.ref, refseq, lastpos, curpos, '-'):
            if is_hpoly(hpolys, col.ref, rec.POS):
                rec.INFO['HPOLY'] = True
            out_vcf.write_record(rec)
        # Generate the vcf row for that column
        row = generate_vcf_row(col, refseq, minbq, maxd, mind, minth, biasth, bias)
        if is_hpoly(hpolys, col.ref, curpos):
            if row.INFO['CB'] == 'N':
                row = generate_vcf_row(col, refseq, 10, maxd, 2, 0.5, biasth, bias)
            row.INFO['HPOLY'] = True
        # Write the record to the vcf file
        out_vcf.write_record(row)
        # Set last position seen
        lastpos = row.POS

    # Insert blank vcf records from last position in mpileup to the end of regionstring
    for rec in blank_vcf_rows(refname, refseq, lastpos, refend+1, '-'):
        if is_hpoly(hpolys, refname, rec.POS):
            rec.INFO['HPOLY'] = True
        out_vcf.write_record(rec)

    # Close the file
    out_vcf.close()

    return output_path

def blank_vcf_rows(refname, refseq, frompos, topos, call='-'):
    '''
    Returns a list of blank vcf rows for all positions that are missing
    between frompos and topos.

    Not inclusive of topos or frompos

    The blank rows should represent a gap in the alignment and be called whatever
    call is set too

    :param str refname: Reference name to set _Record.CHROM
    :param str refseq: Reference sequence to get reference base from
    :param int topos: Current position in the alignment(1 based)
    :param int frompos: Last position seen in the alignment(1 based)
    :param str call: - What to set the CB info field to(Default to:)
    :param bool includeend: Include end base position

    @returns a list of vcf.model._Record objects filled out with the DP,RC,RAQ,PRC,CBD=0 and CB=call
    '''
    #print refname
    #print 'Topos: {0}'.format(topos)
    #print 'formpos: {0}'.format(frompos)
    # Blank record list
    records = []
    # Only do records between frompos and topos
    for i in range(frompos + 1, topos):
        # Reference base at current position
        rec = blank_vcf_row(refname, refseq, i, call)
        records.append(rec)

    return records

def blank_vcf_row(refname, refseq, pos, call='-'):
    '''
    Generates a blank VCF row to insert for depths of 0

    :param str refseq: Bio.seq.seq object representing the reference sequence
    :param str pos: Reference position to get the reference base from(1 indexed)
    :param str call: What to set the CB info field to

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
    record = vcf.model._Record(refname, pos, None, refseq[pos-1], '.', None, None, info, None, None)
    return record

def bias_hq(stats, biasth=50, bias=10):
    '''
    Biases high quality reads in stats so that they are more likely to be selected later on.
    Essentially duplicates quality scores(baseq values) by a factor of bias.
    Given the example::

        stats = { 'A': {'baseq':[40,40,40]}, 'C': {'baseq':[50,50,50,50,50,50,50]} }

    Without any biasing the called base would likely be called an M as there is no majority of C's
    since A is 30% and C is 70%(Assuming 0.8 minth which requires 80% to be called)
    With bias_hq(stats, 50, 2) you would end up with 3 A's and 14 C's or 18% A and 82% C so C
    would be more correctly called for the consensus

    bias will round the length of the created list to the nearest integer so if the length of a
    baseq is say 10 and the bias is 1.25 the resulting baseq list will be 13 in length(addition of 3 baseq)

    :param dict stats: Should most likely be a stats2 dictionary although stats would work as well
    :param int biasth: What quality value(>=) should be considered to be bias towards
    :param int bias: How much to bias aka, how much to multiply the # of quals >= biasth(has to be int >= 1)

    :rtype: dict
    :return: stats2 formatted dictionary with all baseq lists appended to with the bias amount
    '''
    if bias < 1 or int(bias) != bias:
        raise ValueError("bias was set to {} which is less than 1. Cannot bias on a factor < 1".format(bias))

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
        stats2[k]['mapq'] = v.get('mapq', [])
        stats2['depth'] += len(stats2[k]['baseq'])
    return stats2

def pile_stats(mpileupcol, refbase, minbq, mind, biasth, bias):
    '''
    Returns the modified statistics from an mpileupcol that are suitable for base calling

    All parameters are identical to generate_vcf_row

    @returns a stats2 dictionary that is modified by biasing reference bases and high quality bases
    '''
    # This call needs to be replaced after the mpileupcol parameter is put into place
    # then stats_for_col(mpileup
    s = mpileupcol.base_stats()
    # Bias high quality first as it may change the behavior of mark_lq as the depth may
    # increase above the mind threshold
    stats2 = bias_hq(s, biasth, bias)
    stats2 = mark_lq(stats2, minbq, mind, refbase)

    # Update stats2 so that it does not include low quality bases since we
    # are equal to or above the min depth
    if stats2['depth'] >= mind:
        if '?' in stats2:
            stats2['depth'] -= len(stats2['?']['baseq'])
            del stats2['?']

    return stats2

def generate_vcf_row(mpileupcol, refseq, minbq, maxd, mind=10, minth=0.8, biasth=50, bias=10):
    '''
    Generates a vcf row and returns it as a string

    :param str mpileupcol: samtools.MpileupColumn
    :param str refseq: Bio.Seq.Seq object representing the reference sequence
    :param str minbq: minimum base quality to be considered or turned into an N
    :param str maxd: Maximum depth for pileup
    :param str mind: Minimum depth decides if low quality bases are N's or if they are removed
    :param str minth: Minimum percentage to call a base(unless no bases have > minth then the maximum pct base would be called
    :param int biasth: What quality value(>=) should be considered to be bias towards
    :param int bias: How much to bias aka, how much to multiply the # of quals >= biasth(has to be int >= 1)

    @returns a vcf.model._Record
    '''
    from collections import OrderedDict
    # The base position should be the same as the second item in the parsed region string
    start = mpileupcol.pos
    # Get the reference base from the reference sequence
    # Python is 0-index, biology is 1 index
    rb = refseq[start-1].upper()

    stats2 = pile_stats(mpileupcol, rb, minbq, mind, biasth, bias)

    # info needs to contrain the depth, ref count #, % ref count, Ave ref qual, alt count #, % ref count, Ave alf qual
    # Holds the info dictionary in order
    info = {}
    # Alternate base stats info(PAC,AC,AAQ...)
    alt_info = info_stats(stats2, rb)
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
        info['RAQ'] = int(round(sum(refstats['baseq']) / float(len(refstats['baseq'])), 0))
        # Percentage Reference Count is len of qualities / depth
        info['PRC'] = int(round((100.0 * len(refstats['baseq'])) / float(stats2['depth']), 0))
    else:
        info['RC'] = 0
        info['RAQ'] = 0
        info['PRC'] = 0

    # Get Called Base Info
    cb, cbd = caller(stats2, minbq, maxd, mind, minth)
    # Re-evaluate an N call if gaps are involved
    if '*' in stats2 and cb == 'N':
        cb, cbd = caller(stats2, minbq, maxd, mind, 0.51)

    # Set the Called base
    info['CB'] = cb
    info['CBD'] = cbd

    if not alt_bases:
        alt_bases = '.'

    # need to record each line of the vcf file.
    record = vcf.model._Record(mpileupcol.ref, start, None, rb, alt_bases, None, None, info, None, None)

    return record

def caller(stats2, minbq, maxd, mind=10, minth=0.8):
    '''
    Calls a given base at refstr inside of bamfile. At this time refstr has to be a single
    base position('refname':N-N). The base is determined by first labeling all bases less than minbq as N and then
    determining if the depth is < mind or >= mind.
    
    If < and the % of N is > minth then call it an N as it is the majority.
    If >= 10 then remove all N

    The final stage is to call call_on_pct with the remaining statistics

    :param str stats: stats dictionary generated from stats_at_refpos.stats
    :param str minbq: Minimum base quality to determine if it is low quality or not to call it an N
    :param str maxd: Maximum depth in the pileup to allow
    :param str mind: Minimum depth threshold
    :param str minth: Minimum percentage for an base to be called non ambigious

    @returns one of the items from the set('ATGCMRWSYKVHDBN')
    '''
    # Else we will determine if the N's are the majority now
    # defines if the base is an N
    if '?' in stats2:
        nlen = len(stats2['?']['baseq'])
        np = nlen/(stats2['depth']*1.0)
        if np > (1-minth):
            return ('N', nlen)
    return call_on_pct(stats2, minth)

def call_on_pct(stats2, minth=0.8):
    '''
    Calls a base from the given stats dictionary if it is the majority. A majority base is
    any base where it is in %total >= minth.

    Base quality is not used to make the determination, but instead the number of quality scores
    in the baseq list is used as it depicts the depth of that base. The quality scores are not used
    as the stats dictionary should already be run through the label_n function.

    :param str stats2: Stats dictionary returned from mark_lq or stats_at_refpos.stats
    :param str minth: minimum percentage that a base needs to be present in order to be called non-ambiguous

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
            if np_2 > round((1-minth),2):
                nt_list += base
                count += len(bquals)
    dnalist = sorted(nt_list)
    try:
        return (iupac_amb(dnalist), count)
    except ValueError as e:
        if dnalist:
            ndepth = count
        else:
            ndepth = stats2.get('depth', 0)
        return ('N', ndepth)

def info_stats(stats, rb):
    '''
    Returns a dictionary that can be used to fill in the info field of a vcf record
    The returned dictionary will contain the AC(Alternate Count), AAQ(Alternate Average Quality)
    PAC(Percentage Average Count) and bases.

    AC should be a list of integers
    AAQ should be a list of integers
    PAC should be a list of integers rounded to the nearest integer
    bases should be the ordered list of bases for each item in AC,AAQ, PAC such that
    zip(bases, AC/AAQ/PAC) would work as expected

    :param str stats: stats dictionary(should accept either stats or stats2)
    :param str rb: Reference base to ignore in the outputted dictionary

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
