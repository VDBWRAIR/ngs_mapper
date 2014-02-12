from stats_at_refpos import stats
from alphabet import iupac_amb

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

        @param bamfile - 

        @returns path to vcf_output_file
    '''
    pass

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
    # call the nucleotides
    #s = stats(bamfile, regionstr, 1, 1, maxd)
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
