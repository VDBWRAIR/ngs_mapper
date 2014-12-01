========
Analysis
========

Running the Pipeline
====================

New Run Setup
-------------

#. Ensure the run you want to work with has been synced to the storage server
#. Create a :doc:`samplesheet`
    For standardization purposes, the samplesheet name should probably be samplesheet.tsv

Single Sample
-------------

If you need to run a single sample you can now use the :doc:`scripts/runsample.py` script to do so easily.
Things you need to provide :doc:`runsample.py`:

* Path to the ReadsBySample directory for the samplename
    * Example: /path/to/NGSData/ReadsBySample/1090-01
* Path to the reference file
    * Example: /path/to/Analysis/References/Den1__WestPac__1997.fasta

Run the sample
--------------

#. Setup the environment

    .. code-block:: bash

        . $HOME/.miseqpipeline/bin/activate

#. Run the sample(this will put the project in the current directory you are in, you may want to look at the -od argument of :doc:`runsample.py`)

    .. code-block:: bash

        SAMPLE=samplename
        runsample.py /path/to/NGSData/ReadsBySample/${SAMPLE} /path/to/reference ${SAMPLE} -od Projects/${SAMPLE}

    You need to replace samplename with the name of your sample and also replace /path/to/reference with the actual path to the reference you want to use

Changing defaults for pipeline stages
=====================================

If you want to change any of the settings of any of the pipeline stages you will need to create a :doc:`config.yaml` and supply it to :doc:`runsample.py` using the -c option. You can read more about how to create the config and edit it via the :doc:`scripts/make_example_config` script's page

Multiple Samples
================

Typically you will want to run many samples in parallel by utilizing a :doc:`samplesheet`. These instructions show how to do so.

#. Setup the environment

    .. code-block:: bash

        . $HOME/.miseqpipeline/bin/activate

#. Run the samples

    .. code-block:: bash

        runsamplesheet.sh /path/to/NGSData/ReadsBySample samplesheet.tsv

    * **Note** If you did not name the samplesheet samplesheet.tsv you will need to put the path to it instead of ../samplesheet.tsv in the above command
    * Basically just calls :doc:`scripts/runsample.py` over and over for each sample/reference pair in the :doc:`samplesheet` that you created

Rerunning Samples
=================

Rerunning samples is very similar to just running samples.

#. Copy and edit the existing :doc:`samplesheet` and comment out or delete the samples you do not want to rerun.
#. Run the :doc:`scripts/runsamplesheet.sh` script on the modified samplesheet
    * **Note**: As of right now, you will have to manually remove the existing project directories that you want to rerun.
#. Regenerate graphics for all samples
    * The -norecreate tells it not to recreate the qualdepth.json for each sample which is very time consuming. The reran samples should already have recreated their qualdepth.json files when runsample.py was run on them.

        .. code-block:: bash

            graphs.sh -norecreate

#. You should not have to rerun :doc:`scripts/consensuses.sh` as it just symlinks the files
