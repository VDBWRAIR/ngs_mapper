==============
Data Structure
==============

At this time data is organized inside of what is called the NGS Data Structure. This structure is composed of 3 critical directories.

* RawData
* ReadData
* ReadsBySample

Getting Data into the data structure
====================================

See :doc:`ngsdatasync`

Diagram
=======

.. image:: _static/ngsdata_diagram.png

RawData
=======

RawData is composed of all files that originate from each of the instruments' Run.
Some instruments may create ReadData as well or very close to ReadData, but it is still considered RawData.

Some examples of RawData would be:

* Run_3130xl\_ directories containing \*.ab1 files(Sanger)
* Directories under the MiSeqOutput directory(MiSeq)
* R\_\* directories containing signalProcessing or fullProcessing directories(Roche)

ReadData
========

ReadData is any sequence file format that can be utilized by NGS mapping/assembly applications.
At this time these file formats typically end with the following extensions:

    .ab1
    .sff
    .fastq

ReadsBySample
=============

This directory contains only directories that are named after each of the samplenames that have been sequenced. The concept of this folder is to make it very easy to look up all data related to a given samplename.

The files that live in these directories must adhere strictly to a specific naming scheme so the pipeline can determine what reads come from which platforms as they need to be handled slightly different.

See the following naming regular expressions defined in :py:mod:`ngs_mapper.data` for more information about how these filenames have to be named:

* :py:mod:`sanger <ngs_mapper.data.SANGER>`
* :py:mod:`miseq <ngs_mapper.data.MISEQ>`
* :py:mod:`roche <ngs_mapper.data.ROCHE>`
* :py:mod:`iontorrent <ngs_mapper.data.IONTORRENT>`

If you have files that do not match these naming standards the pipeline will essentially ignore them and you may get errors when you run :py:mod:`runsample.py <ngs_mapper.runsample>`. Inspecting the logs you will see errors about Somehow no reads being compiled which is a good indication that your files are named incorrectly.

To ensure your files are named correctly you can use the various sync scripts listed in :doc:`ngsdatasync`
