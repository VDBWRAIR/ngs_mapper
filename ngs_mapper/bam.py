import subprocess
import samtools
import filehandle

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
    p = subprocess.Popen( cmd, stdin=input, stdout=output )

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
    p = subprocess.Popen( cmd, stdin=input )
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

    p = subprocess.Popen( cmd )
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

    p = subprocess.Popen( cmd )
    p.wait()

    return sortedbam + '.bai'

def get_refstats( bamfile ):
    '''
        UNTESTED
        Run samtools idxstats on the given bamfile

        @returns dictionary keyed by refname and values of [refname, reflen, #mapped reads, singletons]
    '''
    cmd = ['samtools','idxstats',bamfile]
    p = subprocess.Popen( cmd, stdout=subprocess.PIPE )
    sout,serr = p.communicate()
    return {line.split()[0]:line.split() for line in sout.splitlines()}

def bam_to_fastq(input_fh):
    '''
    Convert a bam file to fastq format by simply extracting the first, 10th and 11th
    columns from the sam output

    http://samtools.github.io/hts-specs/SAMv1.pdf
    1st column is QNAME
    10th column is SEQ
    11th column is QUAL

    :param str|file input_fh: path to bam or file handle
    :return: generator object yielding fastq string records
    '''
    # Normalize to make sure we end up with an already opened
    # filehandle(gzip files should work too since filehandle.open does that
    if isinstance(input_fh, str):
        fh, ext = filehandle.open(input_fh)
    else:
        fh = input_fh

    # samrows is just an stdout iterator
    samrows = samtools.view(fh)

    # fastq format string
    fqformat = '@{read}\n{seq}\n+\n{qual}'

    # Iterate the rows and format to fastq
    for samrow in samrows:
        row = samrow.split('\t')
        qname = row[0]
        seq = row[9]
        qual = row[10]
        yield fqformat.format(
            read = qname,
            seq = seq,
            qual = qual
        )
