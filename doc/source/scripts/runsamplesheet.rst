=================
runsamplesheet.sh
=================

Runs :py:mod:`runsample <miseqpipeline.runsample>` on every sample/reference pair inside of a :doc:`samplesheet`

Usage
=====

.. code-block:: bash

    runsamplesheet.sh /path/to/ReadsBySample /path/to/samplesheet.tsv

Passing options to runsample.py
-------------------------------

You can run runsamplesheet.sh and pass arguments to runsample.py by prepending RUNSAMPLEOPTIONS="" to the command

Example: adding -minth option
-----------------------------

This would run each sample and pass "-minth 0.95" to :py:mod:`runsample <miseqpipeline.runsample>`

.. code-block:: bash

    RUNSAMPLEOPTIONS="-minth 0.95" runsamplesheet.sh /path/to/ReadsBySample /path/to/samplesheet.tsv

Creates
=======

* Directory named Projects
* All output from :py:mod:`runsample <miseqpipeline.runsample>` placed under Projects
* Pipeline Graphics
    * :py:mod:`sample_coverage <miseqpipeline.coverage>`
    * :doc:`graphs`
