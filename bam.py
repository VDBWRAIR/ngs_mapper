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
