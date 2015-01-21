=======
Scripts
=======

User Scripts
------------

These are scripts that you run directly that will run the Supplemental scripts

* :py:mod:`stats_at_refpos <ngs_mapper.stats_at_refpos>`
* :doc:`runsamplesheet`
* :py:mod:`runsample <ngs_mapper.runsample>`
* :doc:`graphs`
* :doc:`consensuses`
* :py:mod:`miseq_sync <ngs_mapper.miseq_sync>`
* :py:mod:`roche_sync <ngs_mapper.roche_sync>`
* :py:mod:`sanger_sync <ngs_mapper.sanger_sync>`
* :py:mod:`ion_sync <ngs_mapper.ion_sync>`
* :py:mod:`rename_sample <ngs_mapper.rename_sample>`
* :py:mod:`make_example_config <ngs_mapper.config>`

Supplemental
------------

These are scripts that you can run manually, however, they are run automatically by the User Scripts above

* :py:mod:`run_bwa_on_samplename <ngs_mapper.run_bwa>`
* :py:mod:`vcf_consensus <ngs_mapper.vcf_consensus>`
* :doc:`gen_flagstats`
* :py:mod:`graphsample <ngs_mapper.graphsample>`
* :py:mod:`graph_mapunmap <ngs_mapper.graph_mapunmap>`
* :py:mod:`tagreads <ngs_mapper.tagreads>`
* :py:mod:`base_caller <ngs_mapper.base_caller>`
* :py:mod:`graph_times <ngs_mapper.graph_times>`
* :py:mod:`trim_reads <ngs_mapper.trim_reads>`
* :py:mod:`fqstats <ngs_mapper.fqstats>`
* :py:mod:`sample_coverage <ngs_mapper.coverage>`

Libraries
---------

Python Scripts/Modules that you can import to do other analysis

* :py:mod:`ngs_mapper.run_bwa`
* :py:mod:`ngs_mapper.reads`
* :py:mod:`ngs_mapper.data`
* :py:mod:`ngs_mapper.bam`
* :py:mod:`ngs_mapper.alphabet`
* :py:mod:`ngs_mapper.stats_at_refpos`
* :py:mod:`ngs_mapper.samtools`
* :py:mod:`ngs_mapper.log`

Deprecated
----------

Scripts that are no longer used, but kept for reference in the deprecated directory

* varcaller.py
* variants.sh
* perms.sh
* gen_consensus.sh
* setup
* install.sh
