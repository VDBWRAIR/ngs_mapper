from subprocess import Popen, PIPE
import numpy as np
import itertools

def mpileup( bamfile, regionstr=None, minmq=20, minbq=25, maxd=100000 ):
    '''
        A simple  wrapper around the samtools mpileup command to ensure that
        samtools mpileup is called the way we would expect for MPilupColumn to work

        @param bamfile - path to a bam file
        @param regionstr - Region string acceptable to the -r option for mpileup. If None is provided,
            then the same output as not specifying the -r option to mpileup is expected
        @param minmq - Minimum mapping qualty threshold. Same as -q option to mpileup
        @param minbq - Minimum base quality or min BAQ. Same as -Q option to mpileup
        @param maxd - Maximum depth to consider. Same as -d option to mpileup

        @returns the stdout file descriptor handle
    '''
    # The command that will be executed
    cmd = ['samtools','mpileup','-s','-q',str(minmq),'-Q',str(minbq),'-d','{}'.format(maxd)]
    # Only include -r if regionstr is set
    if regionstr:
        cmd += ['-r',regionstr]
    # Add the bamfile to the last parameter
    cmd.append( bamfile )
    # The command will be executed and output sent through the pipe so that it can be iterated
    # over
    p = Popen( cmd, stdout=PIPE, stderr=open('/dev/null','w') )
    # Return the stdout file descriptor handle so it can be easily iterated
    return p.stdout

def char_to_qual( qual_char ):
    '''
        Converts a given quality character to the phred - 33 integer
        
        @param qual_char - Quality character to convert to a phred - 33 integer quality

        @return phred - 33 quality integer
    '''
    return ord( qual_char ) - 33

class MPileupColumn(object):
    '''
        Represents a single Mpileup column

        Object is easily initialized by supplying the constructor with an mpileup string.
        Note that if older samtools is used that does not contain the fix with respect to mapping quality length != base quality length
        then the bquals list will be None

        It is assumed that minimum base quality and mapping quality have already been applied to the mpileup string that is 
        provided to the constructor.

        @param mpileup_str - Mpileup string
    '''
    _bases = ''
    _mquals = ''
    _bquals = ''
    def __init__( self, mpileup_str ):
        parts = mpileup_str.rstrip().split('\t')
        if len(parts) == 7:
            self.ref,self.pos,self.refbase,self.depth,self._bases,self._bquals,self._mquals = parts
        else:
            self.ref,self.pos,self.refbase,self.depth,self._bases,self._bquals = parts
            self._mquals = ''

    @property
    def depth( self ):
        return self.__dict__['depth']
    
    @depth.setter
    def depth( self, value ):
        ''' Depth is an integer '''
        self.__dict__['depth'] = int(value)

    @property
    def pos( self ):
        return self.__dict__['pos']
    
    @pos.setter
    def pos( self, value ):
        ''' Position is an integer '''
        self.__dict__['pos'] = int(value)

    @property
    def bases( self ):
        '''
            Returns the bases with the inserts, deletions, $ and ^qual removed.
            This means it returns just the bases that are really of interest.
            it also includes the * which indicates a deletion.
        '''
        # Will contain the cleaned base
        cleaned = ''
        indel = False
        indel_count = None
        indel_base = ''
        for x in self._bases:
            # If a new insert or delete is starting
            if x in '+-':
                indel = x
                indel_count = None
                continue
            # Only if we are parsing through an indel
            if indel != False:
                # Parse the integer of how many bases to chomp
                if indel in '+-' and indel_count is None:
                    indel_count = int(x)
                    continue
                # We are on a base that needs to be chomped
                if indel in '+-':
                    # Chomp the base
                    if indel_count > 0:
                        indel_count -= 1
                        indel_base = x
                        continue
                    # Last base in the series
                    else:
                        if indel == '+':
                            # Indel bases do not have base quality nor mapping quality
                            # so do not include them in self.bases
                            #cleaned += indel_base
                            pass
                        indel = False
                        indel_count = None
            # Add any normal base, but make it uppercase
            if x.upper() in ('ACTGN*'):
                cleaned += x.upper()
            # . and , mean a match to the reference base
            elif x in ('.,'):
                cleaned += self.refbase
        return cleaned

    @property
    def bquals( self ):
        '''
            Returns the base qualities as a phred - 33 integer
        '''
        return map(char_to_qual, self._bquals)

    @property
    def mquals( self ):
        '''
            Returns the mapping qualities as a phred - 33 integer
            Truncates mquals to the same length as base qualities as older samtools
            had a bug where it would return all mapping qualities regardless of the -q and -Q 
            thresholds being set.
            It will return the empty list if it would otherwise truncate due to the reason above, but
            all the values are not the same since there would be no way to tell what qual values
            match what bases.
        '''
        # Check to make sure map qual len is same as base qual length
        if len(self._bquals) == len(self._mquals):
            return map(char_to_qual, self._mquals)
        # Otherwise we can only proceed if all items are the same
        elif len(set(self._mquals)) == 1:
            l = len(self._bquals)
            return map(char_to_qual, self._mquals[:l])
        else:
            return []

    def bqual_avg( self ):
        ''' Returns the mean of the base qualities rounded to 2 places '''
        return round( np.mean( self.bquals ), 2 )

    def mqual_avg( self ):
        ''' Returns the mean of the mquals rounded to 2 places '''
        return round( np.mean( self.mquals ), 2 )

    def __iter__( self ):
        '''
            Returns an iterator that zips together bases, base qualities and mapping qualities
            Since mquals sometimes may be missing izip_longest will fill it with 0's
        '''
        return itertools.izip_longest( self.bases, self.bquals, self.mquals, fillvalue=0 )

    def base_stats( self ):
        '''
            Returns a compiled statistics dictionary for all the different base values in this column
            The dictionary must contain the following keys:
                depth: Total depth of this column which should be self.depth
                bqualsum: sum of the base qualities
                mqualsum: sum of the mapping qualities
                'A/C/T/G/N/*': dictionary of information about the mapping and base qualities for an individual base
                    mapq: list of all the mapping qualities for this base(int(phred - 33))
                    baseq: list of all the base qualities for this base(int(phred - 333))

            @returns the stats dictionary
        '''
        bquals = self.bquals
        mquals = self.mquals
        bqualsum = float( sum( bquals ) )
        mqualsum = float( sum( mquals ) )
        # Lets just make sure of a few things because samtools mpileup isn't exactly documented the best
        assert len(bquals) == self.depth, "Somehow length of bases != length of Base Qualities"
        depth = self.depth
        stats = {'depth':depth,'mqualsum':mqualsum,'bqualsum':bqualsum}
        for b,bq,mq in self:
            if b not in stats:
                stats[b] = {'baseq':[],'mapq':[]}
            stats[b]['baseq'].append(bq)
            stats[b]['mapq'].append(mq)
        return stats

    def __str__( self ):
        ''' Returns the mpileup string '''
        return "{ref}\t{pos}\t{refbase}\t{depth}\t{_bases}\t{_bquals}\t{_mquals}".format(**self.__dict__)
