from subprocess import Popen, PIPE
import numpy as np
import itertools
import re

def view( infile, *args, **kwargs ):
    '''
        A simple wrapper around samtools view command that will just return the stdout iterator
        The command should accept the same arguments as samtools view except that any single argument
        such as S, h or b just need to be kwargs set to true(S=True,h=True,b=True) otherwise their value
        would just be the same as what the samtools view command accepts

        @param infile - The sam|bam file to read or an open file-like object to read
        
        @returns file like object representing the output of view
    '''
    # The base command
    cmd = ['samtools','view']
    # Put in all the kwargs
    for k,v in kwargs.items():
        # If value is a bool then the k is a flag argument with no value
        if v is True:
            cmd.append( '-'+k )
        else:
            cmd.append( '-'+k )
            cmd.append( v )
    # Now put the file we are working with or - for stdin
    if isinstance(infile,str):
        cmd.append( infile )
        in_handle = False
    else:
        cmd.append('-')
        in_handle = infile

    # Were there any regionstrings specified?
    cmd += args

    # Execute the command
    if in_handle:
        p = Popen( cmd, stdout=PIPE, stdin=in_handle.fileno() )
    else:
        p = Popen( cmd, stdout=PIPE )

    # Return the stdout iterator
    return p.stdout

class Prop(object):
    '''
    Defines a property that auto converts the value
    to the defined type
    '''
    def __init__( self, name, type=int ):
        self.name = name
        self.type = type

    def __get__( self, obj, objtype ):
        return obj.__dict__[self.name]

    def __set__( self, obj, val ):
        obj.__dict__[self.name] = self.type(val)

class SamRow(object):
    '''
    Represents a single sam row
    
    Object is instantiated by supplying it with a valid sam row string

    @param samrow_str - Sam row string
    '''
    FLAG = Prop( 'FLAG', int )
    MAPQ = Prop( 'MAPQ', int )
    POS = Prop( 'POS', int )
    TLEN = Prop( 'TLEN', int )
    PNEXT = Prop( 'PNEXT', int )

    def __init__( self, samrow_str ):
        # Only split up to 11 times. The last element will be all the tags if they are there at all
        parts = samrow_str.rstrip().split( '\t', 11 )
        if len( parts ) == 12:
            self.QNAME, self.FLAG, self.RNAME, self.POS, self.MAPQ, self.CIGAR, \
                self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, self._qual, self._tags = parts
        else:
            self.QNAME, self.FLAG, self.RNAME, self.POS, self.MAPQ, self.CIGAR, \
                self.RNEXT, self.PNEXT, self.TLEN, self.SEQ, self._qual = parts
            self._tags = ''

    @property
    def TAGS( self ):
        '''
        Returns python objects for each flag type
        '''
        tags = re.findall( '([A-Za-z]{2}):([AifZHB]):(\S+)', self._tags )
        t = []
        for name,typ,val in tags:
            if typ == 'A':
                value = str(val)
            elif typ == 'i':
                value = int(val)
            elif typ == 'f':
                value = float(val)
            elif typ == 'Z':
                value = str(val)
            elif typ == 'H':
                value = hex(int(val,0))
            elif typ == 'B':
                value = [int(x) for x in val.split(',')]
            else:
                raise ValueError("{0} is not a supported field type".format(typ))

            t.append( (name,value) )
            
        return t

    @property
    def QUAL( self ):
        '''
        Returns the quality scores as a list of integers
        '''
        return [char_to_qual(c) for c in self._qual]

    def __str__( self ):
        s = '{QNAME}\t{FLAG}\t{RNAME}\t{POS}\t{MAPQ}\t{CIGAR}\t{RNEXT}' \
                '\t{PNEXT}\t{TLEN}\t{SEQ}\t{_qual}'.format(**self.__dict__)
        if self._tags:
            s += '\t'+self._tags
        return s

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

    @returns file like object representing the output of mpileup
    '''
    # The command that will be executed
    cmd = ['samtools','mpileup','-s','-q',str(minmq),'-Q',str(minbq),'-d','{0}'.format(maxd)]
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

def nogap_mpileup(*args, **kwargs):
    '''
    Wrapper around mpileup that fills in missing positions with 0 depth
    Arguments are the same as mplileup

    Returns a generator of mpileup rows
    '''
    lastref = None
    lastpos = 0
    for pile in mpileup(*args, **kwargs):
        col = pile.split('\t')
        refname = col[0]
        pos = int(col[1])
        # First iteration
        if lastref is None:
            lastref = refname
        elif lastref != refname:
            # New reference
            lastref = refname
            lastpos = 0

        # yield all blank rows needed
        for i in range(lastpos+1, pos):
            row = '{0}\t{1}\t\t0\t\t\t'.format(
                refname, i
            )
            yield row

        yield pile
        lastpos = pos

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
        parts = mpileup_str.rstrip('\n').split('\t')
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
        # The counter
        i = 0
        # Iterate every position
        while i < len( self._bases ):
            x = self._bases[i]
            # Just skip - and +
            if x in '-+$':
                i += 1
                continue
            # Skip this position and the next quality scores
            if x == '^':
                i += 2
                continue
            # Add any normal base, but make it uppercase
            if x.upper() in ('ACTGN*'):
                cleaned += x.upper()
                i += 1
            # . and , mean a match to the reference base
            elif x in ('.,'):
                cleaned += self.refbase
                i += 1
            # We are on an indel number now so can just skip them
            else:
                try:
                    # Build the integer that will tell us how many
                    # indel bases to skip
                    n = ''
                    while True:
                        # Just try it if it fails the try except will grab it
                        int(x)
                        n += x
                        # Increment i
                        i += 1
                        # Get new base
                        x = self._bases[i]
                except ValueError as e:
                    i += int( n )
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

            * depth: Total depth of this column which should be self.depth
            * bqualsum: sum of the base qualities
            * mqualsum: sum of the mapping qualities
            * 'A/C/T/G/N/\*': dictionary of information about the mapping and base qualities for an individual base
                * mapq: list of all the mapping qualities for this base(int(phred - 33))
                * baseq: list of all the base qualities for this base(int(phred - 333))

        @returns the stats dictionary
        '''
        bquals = self.bquals
        mquals = self.mquals
        bases = self.bases
        bqualsum = float( sum( bquals ) )
        mqualsum = float( sum( mquals ) )
        # Lets just make sure of a few things because samtools mpileup isn't exactly documented the best
        assert len(bquals) == self.depth, "Somehow length of bases != length of Base Qualities"
        depth = self.depth
        stats = {'depth':depth,'mqualsum':mqualsum,'bqualsum':bqualsum}
        for b,bq,mq in itertools.izip_longest( bases, bquals, mquals, fillvalue=0 ):
            if b not in stats:
                stats[b] = {'baseq':[],'mapq':[]}
            stats[b]['baseq'].append(bq)
            stats[b]['mapq'].append(mq)

        return stats

    def __str__( self ):
        ''' Returns the mpileup string '''
        return "{ref}\t{pos}\t{refbase}\t{depth}\t{_bases}\t{_bquals}\t{_mquals}".format(**self.__dict__)

# Exception for when invalid region strings are given
class InvalidRegionString(Exception): pass

def parse_regionstring( regionstr ):
    '''
        Parses a region string into a 3 item tuple and checks it for errors

        @param regionstr - samtools region string format

        @returns (ref, start, stop)
        @raises InvalidRegionString
    '''
    import re
    if not regionstr:
        raise InvalidRegionString( "{0} is not a valid regionstring".format(regionstr) )
    m = re.match( '(\S+):(\d+)-(\d+)', regionstr )
    if not m:
        raise InvalidRegionString( "{0} is not a valid regionstring".format(regionstr) )

    region = (m.group(1), int(m.group(2)), int(m.group(3)))
    if region[1] > region[2]:
        raise InvalidRegionString( "Start cannot be > stop in a region string: {0}".format(regionstr) )

    return region

