from subprocess import Popen, PIPE

import logging
log = logging.getLogger(__name__)

def samtobam( sam, outbam ):
    '''
        Use samtools to convert a sam file to a bam file
        outbam will be overwritten if it already exists

        @sam - file path or file object(even pipe) of sam file input to convert
        @outbam - file path or file object(even pipe) of bam output destination

        @returns outbam or the file descriptor of the object
    '''
    cmd = ['samtools','view','-Sb','-']
    log.info('Running {}'.format(' '.join(cmd)))
    # Determine if sam is a filepath or file like object
    if isinstance(sam,str):
        input = open(sam)
    else:
        input = sam

    # Determine if outbam is a filepath or file like object
    if isinstance(outbam,str):
        output = open(outbam,'w')
    else:
        output = outbam

    log.debug("CMD: {} STDIN: {} STDOUT: {}".format(cmd,input,output))
    p = Popen( cmd, stdin=input, stdout=output )

    # Returns the original file path or
    # the file like object
    if p.stdout is None:
        # When the output file is an open file object
        # you have to wait for it to finish otherwise nothing will be written
        # when you read it next
        p.wait()
        log.debug("Filename given so returning that filename for output")
        return outbam
    else:
        log.debug("Returning processes stdout value")
        return p.stdout

def sortbam( bam, outbam ):
    '''
        Sorts a bam file using samtools
        outbam cannot be a pipe because samtools index won't allow it
        outbam will be overwritten if it already exists

        @sam - file path or file object(even pipe) of sam file input to convert
        @outbam - file path

        @returns outbam or the file descriptor of the object
    '''
    cmd = ['samtools','sort','-f','-']
    log.info('Running {}'.format(' '.join(cmd)))

    # Determine if sam is a filepath or file like object
    if isinstance(bam,str):
        input = open(bam)
    else:
        input = bam

    # Determine if outbam is a filepath or file like object
    if isinstance(outbam,str):
        cmd.append( outbam )
    else:
        raise ValueError("Output file for sortbam has to be a path not {}".format(outbam))

    log.debug("CMD: {} STDIN: {}".format(cmd,input))
    p = Popen( cmd, stdin=input )
    p.wait()
    return outbam

def mergebams( sortedbams, mergedbam ):
    '''
        Merges the given sortedbams into a file specified by mergedbam

        @param sortedbams - List of sorted bam files to merge(Maybe don't even need to sort them?)
        @param mergedbam - Output file for samtools merge

        @returns the path to mergedbam
    '''
    if not isinstance( sortedbams, list ) or len( sortedbams ) < 2:
        raise ValueError( "Merging bams requires >= 2 bam files to merge. {} was given".format(sortedbams) )

    print sortedbams
    cmd = ['samtools','merge',mergedbam] + sortedbams
    log.info('Running {}'.format(' '.join(cmd)))

    p = Popen( cmd )
    p.wait()

    return mergedbam

def indexbam( sortedbam ):
    '''
        Indexes a sorted bam file

        @sortedbam - file path to a sorted bam file

        @returns the path to the index file for sortedbam(probably sortedbam+'.bai')
    '''
    cmd = ['samtools','index',sortedbam]
    log.info('Running {}'.format(' '.join(cmd)))

    p = Popen( cmd )
    p.wait()

    return sortedbam + '.bai'
