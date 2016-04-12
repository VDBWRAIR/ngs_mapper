====
Help
====

.. toctree::
    :maxdepth: 1
    :hidden:

    Bug Reports <createissue>

Eventually you will run across some errors. No application/software is without bugs.
Here we will compile all of the most common errors and what to look for to find out what is going on

.. _tracebackerror:

Traceback Error
---------------

You will likely encounter a Traceback error at some point due to either a bug or maybe you are running one of the commands incorrectly.

The traceback errors will look like this::

    Traceback (most recent call last):
      File "/home/username/.ngs_mapper/bin/roche_sync", line 9, in <module>
        load_entry_point('ngs_mapper==1.0.0', 'console_scripts', 'roche_sync')()
      File "/home/username/.ngs_mapper/lib/python2.7/site-packages/ngs_mapper/roche_sync.py", line 100, in main
        args = parse_args()
      File "/home/username/.ngs_mapper/lib/python2.7/site-packages/ngs_mapper/roche_sync.py", line 236, in parse_args
        defaults = config['roche_sync']
      File "/home/username/.ngs_mapper/lib/python2.7/site-packages/ngs_mapper/config.py", line 29, in __getitem__
        'Config is missing the key {0}'.format(key)
    ngs_mapper.config.InvalidConfigError: Config is missing the key roche_sync

The easiest way to get good information from the traceback is by working your way backwards(from the bottom to the top).

From this Traceback you should notice that the last line is telling you that the :doc:`config.yaml <config>` file is missing the key roche_sync. You would then edit your config.yaml file and ensure that key exists and then rerun the ``python setup.py install`` portion of the :doc:`install`.

The traceback is simply Python's way of displaying how it got to the error that was encountered. Typically, but not always, the last line of the output contains the most relevant error. If you submit a :doc:`bug report <createissue>`, make sure to include the entire Traceback though.

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

    If you are still not sure, you can search through previous issues on the `GitHub Issue Tracker <https://github.com/VDBWRAIR/ngs_mapper/issues>`_ and/or submit a new :doc:`bug/feature <createissue>`
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
#. There is an error running vcf_consensus that has to do with string index out of bounds
    This has to do with an outdated version of base_caller generating the vcf file you are trying to run vcf_consensus on. See Issue #143 for more information on how to fix that.
#. The pipeline fails on a sample and the log says Somehow no reads were compiled
    This usually indicates that it could not find any reads inside of the location you specified that should contain sample reads. Make sure that the directory you specified when you ran :doc:`scripts/runsamplesheet` or :py:mod:`ngs_mapper.runsample` actually contains a directory with reads for every sample you are running.
    Also check for errors near the top of the log file that say anything about why any reads might have been skipped
#. The pipeline keeps failing on all of my samples or the logs say something about No Space Left On Device
    Please check your /dev/shm and /tmp to see if either is full(``df -h``). You can clear out all of the left-over junk from the pipeline by issuing ``rm -rf /tmp/runsample* /dev/shm/mapbwa*``
    Also, you may need to tell the pipeline to use a different temporary directory. See :ref:`tempdirfiles` for more information.
#. You get a Traceback error that contains ngs_mapper.config.InvalidConfigError: Config is missing the key missingkey
    This indicates that the initial config.yaml file that you created during the :doc:`install` is missing a required key: value pair called missingkey. This most likely happened because you updated the pipeline which introduced new keys in config.yaml.base that you need to add to your config.yaml.
    
    Once you add those new keys, you will need to rerun the ``python setup.py install`` portion of the :doc:`install`.
#. You get errors related to ``no display name and no $DISPLAY environemt variable`` when createing graphics
    See `Issue 75 <https://github.com/VDBWRAIR/ngs_mapper/issues/75>`_
