====
Help
====

.. toctree::
    :maxdepth: 1
    :hidden:

    Bug Reports <createissue>

Eventually you will run across some errors. No application/software is without bugs.
Here we will compile all of the most common errors and what to look for to find out what is going on

.. _faq:

Frequently Asked Questions
--------------------------

#. There is an error. What do I do?
    There are a few log files that you can check. The output on your screen should give you the location of the log file to check for errors.

    As well you can look under the directory of any project and look in files that end in .log

    For instance, if a run fails for any reason it will spit many lines to the screen. When you read through the lines you will see one that mentions "Check the log file" followed by a path to a bwa.log. Navigate to the bwa.log to 
    view a detailed log of what happened.

    There are two other log files which are in the same directory as bwa.log [samplename].std.log and [samplename].log. You can check any of these log files to determine what happend during the run.

    Finally, you can also check the pipeline.log file that is generated when the pipeline is done or if it err'd out.

    If you are still not sure, you can search through previous issues on the `GitHub Issue Tracker <https://github.com/VDBWRAIR/miseqpipeline/issues>`_ and/or submit a new :doc:`bug/feature <createissue>`
#. Where should I run the analysis?
    This is for the most part up to you but eventually you will want the entire analysis folder to end up under /path/to/Analysis somewhere
    You will want to minimize how much the traffic has to travel across the network though. So if you simply create a folder under /path/to/Analysis/PipelineRuns and then you run the pipeline from there, you will essentially be doing the following:

        * Reading the reads across the network for each sample
        * Writing the bam files across the network for each sample
        * Reading the bam files across the network for each sample
        * Writing stats across the network
        * Reading the stats file across the network
        * Writing graphics files across the network

    Suggestion Create the analysis folder somewhere on your computer and run the pipeline there and then transfer the entire folder to the storage server afterwards
#. How many CPUs does my computer have?
    Try running the following command to get how many physical CPU's and how many cores/threads they have

    .. code-block:: bash

        for pid in $(awk '/physical id/ {print $4}' /proc/cpuinfo |sort|uniq)
        do
            echo "--- Processor $pid ---"
            egrep -xA 12 "processor[[:space:]]: $pid" /proc/cpuinfo
        done

#. How many CPUs should I use?
    Check out the command above for more info on how to get how many CPU/Core/Threads you have. Probably best to use (cpu cores \* number of processors)

    If your output was the following then you would probably want to use (2 * 6)

    .. code-block:: bash

        --- Processor 0 ---
        processor         : 0
        ...
        physical id:      : 0
        siblings          : 12
        core id           : 0
        cpu cores         : 6
        ...
        --- Processor 1 ---
        ...
        processor         : 0
        physical id:      : 1
        siblings          : 12
        core id           : 0
        cpu cores         : 6
        ...

    That all being said, you could also try using (number of processors \* siblings) or 24 in the above example,
    but that may actually slow down your analysis
#. How much RAM do I have?
    The following command will tell you how much memory you have in MB

    .. code-block:: bash

        free -m | awk '/Mem:/ {print $2}'

#. The pipeline fails on samples and the bwa.log says something about failing on the reference index
    Make sure to check that you have permissions to read the reference file. The last thing to check is that the reference is formatted correctly in fasta format.
#. There is an error running vcf_consensus.py that has to do with string index out of bounds
    This has to do with an outdated version of base_caller.py generating the vcf file you are trying to run vcf_consensus.py on. See Issue #143 for more information on how to fix that.
#. The pipeline fails on a sample and the log says Somehow no reads were compiled
    This usually indicates that it could not find any reads inside of the location you specified that should contain sample reads. Make sure that the directory you specified when you ran :doc:`scripts/runsamplesheet` or :py:mod:`miseqpipeline.runsample` actually contains a directory with reads for every sample you are running.
    Also check for errors near the top of the log file that say anything about why any reads might have been skipped
#. The pipeline keeps failing on all of my samples or the logs say something about No Space Left On Device
    Please check your /dev/shm to see if it is full as the pipeline uses this special memory filesystem to operate in. You can clear out all of the left-over junk from the pipeline by issuing <pre>rm -rf /dev/shm/runsample* /dev/shm/mapbwa*</pre>
