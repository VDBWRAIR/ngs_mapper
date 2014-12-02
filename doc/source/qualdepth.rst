==============
qualdepth.json
==============

This is the base statistics file that is generated from a given bam file. It is intended to contain all supporting information that will be useful to create graphics and do other statistical analysis through an easy means of simply loading this file.

Typically the qualdepth.json file is generated via the :py:mod:`graphsample.py <miseqpipeline.graphsample>` script, but can be manually created for any given bam file via the BamCoverage's :py:mod:`bam_to_qualdepth.py <miseqpipeline.bam_to_qualdepth>` script.

h2. Loading the file in Python

It is amazingly simple to load json files in python using the build in :py:mod:`json` module

    .. code-block:: python

        import json
        stats = json.load( open('qualdepth.json') )

Now you have access to all of the stats as follows:

* Unmapped Reads

    .. code-block:: python

        stats['unmapped_reads']

* References

    .. code-block:: python

        reflist = [ref for ref in stats if ref != 'unmapped_reads']

* Stats about a given reference

    .. code-block:: python

        ref = reflist[0]
        refstats = stats[ref]
        mapped_reads = refstats['mapped_reads']
        depths = refstats['depths']
        quals = refstats['avgquals']
        reflen = refstats['length']
        minqual = refstats['minq']
        maxqual = refstats['maxq']
        mindepth = refstats['mind']
        maxdepth = refstats['maxd']
