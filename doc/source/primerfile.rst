Primer File
===========

You can set the primer file and the the primer options in the config.yaml file. You'll want to change the fields under trim_reads. To learn more about the options, see the help fields below, and see the documentation at http://www.usadellab.org/cms/?page=trimmomatic.

**Warning:** primer sequences may need to be longer than 20 bases or they don't get trimmed. 

More detailed information can be found in the trimmomatic paper: http://bioinformatics.oxfordjournals.org/content/30/15/2114.full.pdf+html
 
.. code-block:: yaml 

    trim_reads:
        primerfile:
            default:
            help:  'Primer File'
        primerseed:
            default: 2
            help:   'specifies the maximum mismatch count which will still allow a full match to be performed'
        palindromeclip:
            default: 30
            help: 'specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.'
        simpleclip:
            default: 20
            help: 'specifies how accurate the match between any adapter etc. sequence must be against a read.' 

The available options are:

.. code-block:: yaml 

    primerfile: path to the primer file
    seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed
    palindromeClipThreshold: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.
    simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.
