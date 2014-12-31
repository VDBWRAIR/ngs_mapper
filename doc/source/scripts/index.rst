=======
Scripts
=======

User Scripts
------------

These are scripts that you run directly that will run the Supplemental scripts

* :py:mod:`stats_at_refpos <miseqpipeline.stats_at_refpos>`
* :doc:`runsamplesheet`
* :py:mod:`runsample <miseqpipeline.runsample>`
* :doc:`graphs`
* :doc:`consensuses`
* :py:mod:`miseq_sync <miseqpipeline.miseq_sync>`
* :py:mod:`roche_sync <miseqpipeline.roche_sync>`
* :py:mod:`sanger_sync <miseqpipeline.sanger_sync>`
* :py:mod:`rename_sample <miseqpipeline.rename_sample>`
* :py:mod:`make_example_config <miseqpipeline.config>`

Supplemental
------------

These are scripts that you can run manually, however, they are run automatically by the User Scripts above

* :py:mod:`run_bwa_on_samplename <miseqpipeline.run_bwa>`
* :py:mod:`vcf_consensus <miseqpipeline.vcf_consensus>`
* :doc:`gen_flagstats`
* :py:mod:`graphsample <miseqpipeline.graphsample>`
* :py:mod:`graph_mapunmap <miseqpipeline.graph_mapunmap>`
* :py:mod:`tagreads <miseqpipeline.tagreads>`
* :py:mod:`base_caller <miseqpipeline.base_caller>`
* :py:mod:`graph_times <miseqpipeline.graph_times>`
* :py:mod:`trim_reads <miseqpipeline.trim_reads>`
* :py:mod:`fqstats <miseqpipeline.fqstats>`
* :py:mod:`sample_coverage <miseqpipeline.coverage>`

Libraries
---------

Python Scripts/Modules that you can import to do other analysis

* :py:mod:`miseqpipeline.run_bwa`
* :py:mod:`miseqpipeline.reads`
* :py:mod:`miseqpipeline.data`
* :py:mod:`miseqpipeline.bam`
* :py:mod:`miseqpipeline.alphabet`
* :py:mod:`miseqpipeline.stats_at_refpos`
* :py:mod:`miseqpipeline.samtools`
* :py:mod:`miseqpipeline.log`

Deprecated
----------

Scripts that are no longer used, but kept for reference in the deprecated directory

* varcaller.py
* variants.sh
* perms.sh
* gen_consensus.sh
* setup
* install.sh
