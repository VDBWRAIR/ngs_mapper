from subprocess import Popen, PIPE
import numpy as np
import itertools

def mpileup( bamfile, regionstr=None, minmq=20, minbq=25, maxd=100000 ):
    cmd = ['samtools','mpileup','-s','-q',str(minmq),'-Q',str(minbq),'-d','{}'.format(maxd)]
    if regionstr:
        cmd += ['-r',regionstr]
    cmd.append( bamfile )
    p = Popen( cmd, stdout=PIPE )
    return p.stdout

def char_to_qual( qual_char ):
    '''
        Converts a given quality character to the phred - 33 integer
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
        self.__dict__['depth'] = int(value)

    @property
    def pos( self ):
        return self.__dict__['pos']
    
    @pos.setter
    def pos( self, value ):
        self.__dict__['pos'] = int(value)

    @property
    def bases( self ):
        cleaned = ''
        indel = False
        indel_count = None
        indel_base = ''
        for x in self._bases:
            if x in '+-':                                                                                                                                                                                                                                                                                                            
                indel = x
                indel_count = None
                continue
            if indel != False:
                if indel in '+-' and indel_count is None:
                    indel_count = int(x)
                    continue
                if indel in '+-':
                    if indel_count > 0:
                        indel_count -= 1
                        indel_base = x
                        continue
                    else:
                        if indel == '+':
                            #cleaned += indel_base
                            pass
                        indel = False
                        indel_count = None
            if x.upper() in ('ACTGN*'):
                cleaned += x.upper()
            elif x in ('.,'):
                cleaned += self.refbase
        return cleaned

    @property
    def bquals( self ):
        return map(char_to_qual, self._bquals)

    @property
    def mquals( self ):
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
        return round( np.mean( self.bquals ), 2 )

    def mqual_avg( self ):
        return round( np.mean( self.mquals ), 2 )

    def __iter__( self ):
        return itertools.izip_longest( self.bases, self.bquals, self.mquals, fillvalue=0 )

    def base_stats( self ):
        bquals = self.bquals
        mquals = self.mquals
        bqualsum = float( sum( bquals ) )
        mqualsum = float( sum( mquals ) )
        assert len(bquals) == len(mquals), "Somehow length of Base Qualities != Map Qualities"
        depth = self.depth
        stats = {'depth':depth,'mqualsum':mqualsum,'bqualsum':bqualsum}
        for b,bq,mq in self:
            if b not in stats:
                stats[b] = {'baseq':[],'mapq':[]}
            stats[b]['baseq'].append(bq)
            stats[b]['mapq'].append(mq)
        return stats

    def __str__( self ):
        return "{ref}\t{pos}\t{refbase}\t{depth}\t{_bases}\t{_bquals}\t{_mquals}".format(**self.__dict__)
