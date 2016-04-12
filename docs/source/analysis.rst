========
Analysis
========

.. Contents::
    :Depth: 4

Complete Examples
=================

Here we will show you a complete example of running the pipeline using some test data that is included with the source code.

**Note**: Any time you see

.. code-block:: bash

    $> command

It means you should be able to type that into your terminal

All examples will assume your current working directory is inside of the git cloned ngs_mapper directory, aka the following command
ends with ngs_mapper:

.. code-block:: bash

    $> pwd

For both examples below, as always when running the pipeline, you need to ensure 
your installation is in your environments PATH.

See the :doc:`install` for how to setup your PATH


The location of our data sets are under ngs_mapper/tests/fixtures/functional

.. code-block:: bash

    $> ls ngs_mapper/tests/fixtures/functional
    780  780.conf  780.ref.fasta  947  947.conf  947.ref.fasta

Here you can see we have 2 data sets to play with.

* 780 is an H3N2 data set
* 947 is a Dengue 4 data set

You will notice that there is a 780 and a/947 directory
There is also a 780.ref.fasta and 947.ref.fasta file.
The 780 and 947 directory contain all the read files for the 780 and 947 samples
while the 780.ref.fasta and 947.ref.fasta is the reference to map to for each project.
You can ignore the .conf files, they are used by the automated tests.

Quick note about platform identification
----------------------------------------

Reads are identified by the way the first identifier in each file is named.

You can read more about this :ref:`here <platformidentification>`

Using runsample to run a single sample
-----------------------------------------

Some times you just need to run a single sample. Here we will use the `runsample` script to run the 947 example data set and have the analysis be put into a directory called 947 in the current directory.

First, let's see what options there are available for the :py:mod:`runsample <ngs_mapper.runsample>` script to use

This is just an example output and may not match your exact output.

.. code-block:: bash

    $> runsample 
    usage: runsample [-h] [--config CONFIG] [-trim_qual TRIM_QUAL]
                        ...
                        readsdir reference prefix
    runsample: error: too few arguments

What you can take from this is:

* Anything inside of a [] block means that argument to the script is optional and has a default value that will be used if you do not specify it.
* readsdir, reference and prefix are all required arguments that you **MUST** specify

Simplest form of runsample
^^^^^^^^^^^^^^^^^^^^^^^^^^

To run the project with the fewest amount of arguments would be as follows(don't run this, just an example):

.. code-block:: bash

    $> runsample ngs_mapper/tests/fixtures/functional/947 ngs_mapper/tests/fixtures/functional/947.ref.fasta -od 947 947

This will run the 947 data and use the 947.ref.fasta file to map to. All files will be prefixed with 947.
Since we did not specify the -od argument, all the files from the pipeline get dumped into your current directory.

Most likely you will want to specify a separate directory to put all the 947 specific analysis files into. But how?

Getting extended help for runsample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can get extended help information which should print the defualts as well from any script by using the ``--help`` option

.. code-block:: bash

    $> runsample --help
    runsample --help
    usage: runsample [-h] [--config CONFIG] [-trim_qual TRIM_QUAL]
                        [-head_crop HEAD_CROP] [-minth MINTH] [--CN CN]
                        [-od OUTDIR]
                        readsdir reference prefix

    Runs a single sample through the pipeline

    positional arguments:
      readsdir              Directory that contains reads to be mapped
      reference             The path to the reference to map to
      prefix                The prefix to put before every output file generated.
                            Probably the samplename

    optional arguments:
      -h, --help            show this help message and exit
      --config CONFIG, -c CONFIG
                            Path to config.yaml file
      -trim_qual TRIM_QUAL  Quality threshold to trim[Default: 20]
      -head_crop HEAD_CROP  How many bases to crop off the beginning of the reads
                            after quality trimming[Default: 0]
      -minth MINTH          Minimum fraction of all remaining bases after
                            trimming/N calling that will trigger a base to be
                            called[Default: 0.8]
      --CN CN               Sets the CN tag inside of each read group to the value
                            specified.[Default: None]
      -od OUTDIR, --outdir OUTDIR
                            The output directory for all files to be put[Default:
                            /home/myusername/ngs_mapper]

You can see that ``--help`` gives us the same initial output as just running runsample without any arguments, but also contains extended help for all the arguments. The ``--help`` argument is available for all ngs_mapper scripts.
If you find one that doesn't, head over to :doc:`createissue` and file a new Bug Report.

So you can see the -od option's default is our current directory. So if we want our analysis files to go into a specific directory for each sample we run we can specify a different directory. While we are at it, lets try specifying some of the other optional arguments too.

Specifying output directory for analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's tell runsample to put our analysis into a directory called 947 and also tell it to crop off 20 bases from the beginning of each read.

.. code-block:: bash

    $> runsample -od 947 -head_crop 20 ngs_mapper/tests/fixtures/functional/947 ngs_mapper/tests/fixtures/functional/947.ref.fasta 947
    2014-12-22 10:17:52,465 -- INFO -- runsample       --- Starting 947 --- 
    2014-12-22 10:21:28,526 -- INFO -- runsample       --- Finished 947 ---

You can see from the output that the sample started and finished. If there were errors, they would show up in between those two lines and you would have to view the :doc:`help` documentation.

Specifying specific platforms to map
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes you may find the need to only run specific platforms. Maybe you only
will want to run MiSeq read files through the pipeline.

The 947 example project has Roche454, MiSeq and Sanger read files in it, so we
can use it in this example to only map the MiSeq read files

#. Generate your example config which we will edit

    .. code-block:: bash

        make_example_config

#. Now edit the config.yaml file generated in the current directory

   Find the trim_reads section and change the default under platforms to be

        .. code-block:: text

            trim_reads:
                headcrop:
                    default: 0
                    help: 'How many bases to crop off the beginning of the reads after quality
                        trimming[Default: %(default)s]'
                outputdir:
                    default: trimmed_reads
                    help: 'Where to output the resulting files[Default: %(default)s]'
                q:
                    default: 20
                    help: 'Quality threshold to trim[Default: %(default)s]'
                platforms:
                    choices:
                    - MiSeq
                    - Sanger
                    - Roche454
                    - IonTorrent
                    default:
                    - MiSeq
                    #- Sanger
                    #- Roche454
                    #- IonTorrent
                    help: 'List of platforms to include data for[Default: %(default)s]'

    Notice that we have commented out(put # before them) Sanger, Roche454 and IonTorrent.
    You can either comment them out or completely delete them. It is up to you.
#. Then you can run ``runsample`` with the ``-c config.yaml`` argument and it
   will only use MiSeq reads

    .. code-block:: bash

        $> runsample -od 947 -head_crop 20 ngs_mapper/tests/fixtures/functional/947 ngs_mapper/tests/fixtures/functional/947.ref.fasta 947 -c config.yaml

Output from runsample explained
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

So what analysis files were created? You can see them by listing the output directory:

.. code-block:: bash

    $> ls 947
    -rw-r--r--. 1 myusername users 36758279 Dec 22 10:19 947.bam
    -rw-r--r--. 1 myusername users       96 Dec 22 10:19 947.bam.bai
    -rw-r--r--. 1 myusername users    10869 Dec 22 10:21 947.bam.consensus.fasta
    -rw-r--r--. 1 myusername users   269058 Dec 22 10:21 947.bam.qualdepth.json
    -rw-r--r--. 1 myusername users   204502 Dec 22 10:21 947.bam.qualdepth.png
    -rw-r--r--. 1 myusername users  1291367 Dec 22 10:20 947.bam.vcf
    -rw-r--r--. 1 myusername users     2414 Dec 22 10:21 947.log
    -rw-r--r--. 1 myusername users   307180 Dec 22 10:21 947.reads.png
    -rw-r--r--. 1 myusername users    10840 Dec 22 10:17 947.ref.fasta
    -rw-r--r--. 1 myusername users       10 Dec 22 10:18 947.ref.fasta.amb
    -rw-r--r--. 1 myusername users       67 Dec 22 10:18 947.ref.fasta.ann
    -rw-r--r--. 1 myusername users    10744 Dec 22 10:18 947.ref.fasta.bwt
    -rw-r--r--. 1 myusername users     2664 Dec 22 10:18 947.ref.fasta.pac
    -rw-r--r--. 1 myusername users     5376 Dec 22 10:18 947.ref.fasta.sa
    -rw-r--r--. 1 myusername users     2770 Dec 22 10:21 947.std.log
    -rw-r--r--. 1 myusername users    17219 Dec 22 10:18 bwa.log
    -rw-r--r--. 1 myusername users      380 Dec 22 10:20 flagstats.txt
    -rw-r--r--. 1 myusername users      249 Dec 22 10:21 graphsample.log
    -rw-r--r--. 1 myusername users   137212 Dec 22 10:19 pipeline.log
    drwxr-xr-x. 2 myusername users     4096 Dec 22 10:21 qualdepth
    drwxr-xr-x. 2 myusername users     4096 Dec 22 10:18 trimmed_reads
    drwxr-xr-x. 2 myusername users     4096 Dec 22 10:17 trim_stats

You can view information about each of the output files via the :ref:`runsample-output-directory`

Viewing bam files
^^^^^^^^^^^^^^^^^

An easy way to view your bam file quickly from the command line if you have `igv <http://www.broadinstitute.org/igv/>`_  installed is like this:

.. code-block:: bash

    igv.sh -g 947/947.ref.fasta 947/947.bam

Using runsamplesheet.sh to run multiple samples in parallel
-----------------------------------------------------------

:doc:`scripts/runsamplesheet` is just a wrapper script that makes running :py:mod:`runsample <ngs_mapper.runsample>` on a bunch of samples easier.

You just have to first create a :doc:`samplesheet` then you just have to run it as follows:

.. code-block:: bash

    $> runsamplesheet.sh /path/to/NGSData/ReadsBySample samplesheet.tsv

So let's run the 947 and 780 samples as our example.

#. Make a directory for all of our analysis to go into

    .. code-block:: bash

        $> mkdir -p tutorial
        $> cd tutorial

#. Create a new file called samplesheet.tsv and put the following in it(you can use ``gedit samplesheet.tsv`` to edit/save the file)::

    947 ../ngs_mapper/tests/fixtures/functional/947.ref.fasta
    780 ../ngs_mapper/tests/fixtures/functional/780.ref.fasta

#. Run your samplesheet with runsamplesheet.sh

    .. code-block:: bash

        $> runsamplesheet.sh ../ngs_mapper/tests/fixtures/functional samplesheet.tsv
        2014-12-22 12:30:25,381 -- INFO -- runsample       --- Starting 780 --- 
        2014-12-22 12:30:25,381 -- INFO -- runsample       --- Starting 947 --- 
        2014-12-22 12:30:50,834 -- INFO -- runsample       --- Finished 780 ---
        2014-12-22 12:34:08,523 -- INFO -- runsample       --- Finished 947 ---
        1.82user 0.05system 0:01.01elapsed 185%CPU (0avgtext+0avgdata 242912maxresident)k
        0inputs+728outputs (1major+26371minor)pagefaults 0swaps
        5.02user 0.11system 0:04.03elapsed 127%CPU (0avgtext+0avgdata 981104maxresident)k
        0inputs+3160outputs (1major+77772minor)pagefaults 0swaps
        2014-12-22 12:34:19,843 -- WARNING -- graph_times     Projects/780 ran in only 25 seconds
        2014-12-22 12:34:19,843 -- INFO -- graph_times     Plotting all projects inside of Projects

You can see that the pipeline ran both of our samples at the same time in parallel. The pipeline tries to determine how many CPU cores your system has and will run that many samples in parallel.

You can then view all of the resulting output files/directories created

.. code-block:: bash

    $> ls -l
    total 1184
    -rw-r--r--. 1 myusername users   2101 Dec 22 12:34 graphsample.log
    -rw-r--r--. 1 myusername users  50794 Dec 22 12:34 MapUnmapReads.png
    -rw-r--r--. 1 myusername users 756139 Dec 22 12:34 pipeline.log
    -rw-r--r--. 1 myusername users  34857 Dec 22 12:34 PipelineTimes.png
    drwxr-xr-x. 4 myusername users   4096 Dec 22 12:34 Projects
    -rw-r--r--. 1 myusername users 292764 Dec 22 12:34 QualDepth.pdf
    -rw-r--r--. 1 myusername users  52064 Dec 22 12:34 SampleCoverage.png
    -rw-r--r--. 1 myusername users    122 Dec 22 12:28 samplesheet.tsv
    drwxr-xr-x. 2 myusername users   4096 Dec 22 12:34 vcf_consensus

You can view advanced usage and what each of these output files mean by heading over to the :doc:`scripts/runsamplesheet`

Changing defaults for pipeline stages
=====================================

If you want to change any of the settings of any of the pipeline stages you will need to create a :doc:`config` and supply it to :py:mod:`runsample <ngs_mapper.runsample>` using the -c option. You can read more about how to create the config and edit it via the :doc:`config` script's page

Rerunning Samples
=================

Rerunning samples is very similar to just running samples.

#. Copy and edit the existing :doc:`samplesheet` and comment out or delete the samples you do not want to rerun.
#. Run the :doc:`scripts/runsamplesheet` script on the modified samplesheet
    * **Note**: As of right now, you will have to manually remove the existing project directories that you want to rerun.
#. Regenerate graphics for all samples
    * The -norecreate tells it not to recreate the qualdepth.json for each sample which is very time consuming. The reran samples should already have recreated their qualdepth.json files when :py:mod:`runsample <ngs_mapper.runsample>` was run on them.

        .. code-block:: bash

            graphs.sh -norecreate

#. You should not have to rerun :doc:`scripts/consensuses` as it just symlinks the files

.. _tempdirfiles:

Temporary Directories/Files
===========================

The pipeline initially creates a temporary analysis directory for each sample that you run with :py:mod:`runsample <ngs_mapper.runsample>`.

The name of this temporary directory will be samplenameRANDOMrunsample

This directory will be located inside of each project's specified output directory
that was given with ``-od``

If the project fails to complete for some reason then you will need to look inside of that directory for relevant log files to inspect what happened.

Integration with the PBS Schedulers
===================================

runsample has the ability to output a PBS job file instead of running. This may be 
useful if you have access to a PBS Cluster. By default the PBS job that is generated 
is very simplistic.

* The job will change directory to the same directory that qsub is run from
* runsample is then run with the same arguments that were given to generate the
  pbs job without the --qsub arguments.

Example
-------

.. code-block:: bash

    $> runsample ngs_mapper/tests/fixtures/functional/947{,.ref.fasta} 947 --outdir 947test --qsub_l nodes=1:ppn=1 --qsub_M me@example.com
    #!/bin/bash
    #PBS -N 947-ngs_mapper
    #PBS -j oe
    #PBS -l nodes=1:ppn=1
    #PBS -m abe
    #PBS -M me@example.com
    cd $PBS_O_WORKDIR
    runsample ngs_mapper/tests/fixtures/functional/947 ngs_mapper/tests/fixtures/functional/947.ref.fasta 947 --outdir 947test

You can see that the job that was generated essentialy just stripped off any 
--qsub\_ arguments and will rerun the same runsample command in the job.

