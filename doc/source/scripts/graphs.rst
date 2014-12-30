:orphan:

=========
graphs.sh
=========

Purpose
=======

Very simple shell script to create graphics for samples that were mapped.

Creates
-------

* In each sample's directory
    * <samplename>.bam.qualdepth.json
        * Base statistics file for that sample. Read more about [[qualdepth.json]]
    * <samplename>.bam.qualdepth.png
        * Shows Depth/Quality and Mapped/Unmapped for that sample for every reference(Basically all \*.png in qualdepth directory compiled into this file)
    * qualdepth/
        * Directory that contains a qualdepth image for each reference that had reads that mapped
* QualDepth.pdf
    * All qualdepth.png files are compiled into this png to easily view each Sample in one file
    * The QualDepth.pdf file will show you the quality(green) vs. depth(blue) over the mapped genome as well as how many mapped(green) vs. unmapped(red). Typically you want around 37 for the quality and the depth will vary depending on experiment, but you don't want it to dip below 10.
    * You can use the evince command to open the file directly from your terminal ``evince QualDepth.pdf``
* MapUnmapReads.png
    * Shows Mapped/Unmapped reads for each sample as well as Total Mapped/Total Unmapped in one graph to point out samples that have issues
    * This graphic is to show quickly any samples that may have a lot of unmapped reads which could indicate an incorrect reference or some other issues.
* SampleCoverage.png
    * Shows an easy to digest coverage graphic for each sample to know where 'issue' areas are on the genome
    * :py:mod:`sample_coverage <miseqpipeline.coverage>` has more info
* PipelineTimes.png
    * Graphic that shows how many seconds each sample took to run

Basic Usage
===========

Generate all qualdepth.png files inside of each project

    .. code-block:: bash

        graphs.sh

If qualdepth.json files already exist(you already generated them someplace else), then you can run graph.sh and tell it not to recreate them
Very useful if you are simply updating the graphics or recreating them because you accidentally deleted some. If the bam files have changed you will not
want to use this, but instead you will want to recreate the qualdepth.json files(maybe manually so you don't have to recreate all of them which takes a while)

    .. code-block:: bash

        graphs.sh -norecreate
