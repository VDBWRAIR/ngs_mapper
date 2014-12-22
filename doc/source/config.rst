===========
config.yaml
===========

When you install the pipeline you are instructed to copy `miseqpipeline/config.yaml.default <../../../miseqpipeline/config.yaml.default>`_ to miseqpipeline/config.yaml
This file contains all settings that the pipeline will use by default if you do not change them using any of the script options that are available.

When you install the pipeline the config.yaml file gets installed with the pipeline into the installation directory(probably ~/.miseqpipeline)
In order to change the defaults after that you have two options:

* Edit config.yaml inside of the source directory you cloned with git, then go into your miseqpipeline directory and rerun the setup.py command
 
 .. code-block:: bash

        python setup.py install
* Use the :py:mod:`make_example_config <miseqpipeline.config>` to extract the config.yaml into the current directory and use it

Example changing single script defaults
---------------------------------------

If you want to change the quality threshold to use to trim reads when you run :doc:`scripts/trim_reads` you would probably do something as follows:

#. First what options are available for the command?

    .. code-block:: bash

        #> trim_reads.py --help
        usage: trim_reads.py [-h] [--config CONFIG] [-q Q] [--head-crop HEADCROP]
                             [-o OUTPUTDIR]
                             readsdir

        Trims reads

        positional arguments:
          readsdir              Read or directory of read files

        optional arguments:
          -h, --help            show this help message and exit
          --config CONFIG, -c CONFIG
                                Path to config.yaml file
          -q Q                  Quality threshold to trim[Default: 20]
          --head-crop HEADCROP  How many bases to crop off the beginning of the reads
                                after quality trimming[Default: 0]
          -o OUTPUTDIR          Where to output the resulting files[Default:
                                trimmed_reads]

    *You can see that there is a -q option you can specify the quality threshold with*
#. Now run the command with your specific value

    .. code-block:: bash

        #> trim_reads.py -q 5 /path/to/my/input.fastq

This process works pretty slick until you notice that there is no way to easily tell :doc:`scripts/runsample` to specify that same value.
With the version 1.0 release of the pipeline there is now a config file that you can edit and change the Default value any script will use.

Example running :doc:`scripts/runsample` using config.yaml
----------------------------------------------------------

#. First we need to get a config file to work with

    .. code-block:: bash

        #> make_example_config
        /current/working/directory/config.yaml

#. We just need to edit that config.yaml file which should be in the current directory and change the trim_reads's q option default value to 5 then save the file
#. Now just run :doc:`scripts/runsample` as follows

    .. code-block:: bash

        #> runsample.py /path/to/NGSData /path/to/reference.fasta mysample -od mysample -c config.yaml
        2014-11-28 14:39:14,906 -- INFO -- runsample       --- Starting mysample --- 
        2014-11-28 14:39:14,906 -- INFO -- runsample       --- Using custom config from config.yaml ---
        2014-11-28 14:39:35,926 -- INFO -- runsample       --- Finished mysample ---

Example running runsamplesheet.sh using a custom config.yaml
------------------------------------------------------------

You will probably want to be able to run an entire samplesheet with a custom config file as well. If you check out the :doc:`scripts/runsamplesheet` page you will notice that you can specify options to pass on to :doc:`scripts/runsample` by using the RUNSAMPLEOPTIONS variable

#. Generate your config.yaml template

    .. code-block:: bash

        make_example_config

#. Then run :doc:`scripts/runsamplesheet` with your custom config.yaml

    .. code-block:: bash

        #> RUNSAMPLESHEET="-c config.yaml" runsamplesheet.sh /path/to/NGSData/ReadsBySample samplesheet.tsv

Editing config.yaml
===================

The config.yaml file is just a `yaml <http://www.yaml.org>`_ formatted file that is parsed using the python package `pyaml <http://pyyaml.org/>`_
Yaml syntax links for reference:

* `Quick start <http://docs.ansible.com/YAMLSyntax.html>`_
* `More in depth <http://en.wikipedia.org/wiki/YAML>`_

For the miseqpipeline the most important thing is that the NGSDATA value is filled out and contains a correct path to the root of your :doc:`ngsdata`
The rest of the values are pre-filled with defaults that work for most general cases.

Structure of the config.yaml file
---------------------------------

The config.yaml basically is divided into sections that represent defaults for each stage/script that the pipeline has.
It also contains some global variables such as the NGSDATA variable.

Each script/stage requires at a minimum of the default and help defined.

* default defines the default value that option will use
* help defines the help message that will be displayed for that option and probably does not need to be modified
    While yaml does not require you to put text in quotes, it is highly recommended as it will remove some parsing problems if you have special characters in your text such as a : or %
