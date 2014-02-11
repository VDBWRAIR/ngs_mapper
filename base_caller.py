from stats_at_refpos import stats

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



def caller( bamfile, regionstr, minbq, maxd, mind=10, minth=0.8 ):
    '''
        Calls a given base at refstr inside of bamfile. At this time refstr has to be a single
        base position('refname':N-N). The base is determined by first labeling all bases less than minbq as N and then
        determining if the depth is < mind or >= mind.
        If < and the % of N is > minth then call it an N as it is the majority.
        If >= 10 then remove all N

        The final stage is to call call_on_pct with the remaining statistics

        @param bamfile - Path to a bamfile
        @param refstr - Region string acceptable for samtools
        @param minbq - Minimum base quality to determine if it is low quality or not to call it an N
        @param maxd - Maximum depth in the pileup to allow
        @param mind - Minimum depth threshold
        @param minth - Minimum percentage for an base to be called non ambigious

        @returns one of the items from the set( 'ATGCMRWSYKVHDBN' )
    '''
    # call the nucleotides

    s = stats(bamfile, regionstr, 1, 1, maxd)
    stats2 = label_N(s, minbq)


    # if the depth is >= min depth then just remove the N's
    if stats2['depth'] >= mind:
        if 'N' in stats2:        
            del stats2['N']
    else:
        # Else we will determine if the N's are the majority now
        # defines if the base is an N
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

            for q in bquals:
               

                np_2 = len(stats[base]['baseq'])/(stats['depth']*1.0)
                if np_2 > (1-minth):
                    return base


    dnalist = ''.join(sorted(nt_list))


def iupac_amb( dnalist ):
    return (iupac.get(dnalist))
