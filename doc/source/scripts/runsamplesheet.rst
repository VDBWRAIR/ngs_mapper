:orphan:

=================
runsamplesheet.sh
=================

Runs :py:mod:`runsample <ngs_mapper.runsample>` on every sample/reference pair inside of a :doc:`../samplesheet`

Usage
=====

.. code-block:: bash

    runsamplesheet.sh /path/to/ReadsBySample /path/to/samplesheet.tsv

Passing options to runsample
-------------------------------

You can run runsamplesheet.sh and pass arguments to runsample by prepending RUNSAMPLEOPTIONS="" to the command

Example: adding -minth option
-----------------------------

This would run each sample and pass "-minth 0.95" to :py:mod:`runsample <ngs_mapper.runsample>`

.. code-block:: bash

    RUNSAMPLEOPTIONS="-minth 0.95" runsamplesheet.sh /path/to/ReadsBySample /path/to/samplesheet.tsv

Example: Supplying custom config.yaml file
------------------------------------------

#. Generate your custom config.yaml

    .. code-block:: bash

        make_example_config

#. Edit the config.yaml generated to suit your needs
#. Run ``runsamplesheet.sh`` with custom config.yaml

    .. code-block:: bash

        RUNSAMPLEOPTIONS="-c config.yaml" runsamplesheet.sh /path/to/ReadsBySample /path/to/samplesheet.tsv

Creates
=======

* graphsample.log
    * Logfile from running :py:mod:`graphsample <ngs_mapper.graphsample>` on all samples in samplesheet
* MapUnmapReads.png
    * Graphic that shows each sample's mapped vs unmapped read counts
* pipeline.log
    * Logfile that contains essentially the same information on the console you get when you run runsample except it also includes debug lines
* PipelineTimes.png(See :doc:`graphs`)
* Projects
    * All output from :py:mod:`runsample <ngs_mapper.runsample>` placed under Projects named after each sample
* QualDepth.pdf(See :doc:`graphs`)
* SampleCoverage.png(See :doc:`graphs`)
* vcf_consensus
    * Contains symbolic links(shortcuts) to each sample's consensus.fasta file 
