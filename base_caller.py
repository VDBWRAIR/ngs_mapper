def label_N( stats, minbq ):
    '''
        Labels all qualities < minbq as N

        @param stats - Stats dictionary returned from stats_at_refpos.stats
        @param minbq - The mininum base quality to determine if a quality should belong to N

        @returns stats dictionary with baseq subkey for each base in the dictionary.
    '''
    pass

def caller( bamfile, refstr, minbq, maxd, mind=10, minth=0.8 ):
    '''
        Calls a given base at refstr inside of bamfile. At this time refstr has to be a single
        base position. The base is determined by first labeling all bases less than minbq as N and then
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
    pass

def call_on_pct( stats, minth=0.8 ):
    '''
        Calls a base from the given stats dictionary if it is the majority. A majority base is
        any base where it is in %total >= minth.

        @param stats2 - Stats dictionary returned from label_N or stats_at_refpos.stats
        @param minth - minimum percentage that a base needs to be present in order to be called non-ambiguous

        @returns the called base based on the percentages in the given stats
    '''
    pass
