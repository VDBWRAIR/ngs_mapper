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

* Run_3130xl\_ directories containing \*.ab1 files
* Directories under the MiSeqOutput directory
* R\_\* directories containing signalProcessing or fullProcessing directories

ReadData
========

ReadData is any sequence file format that can be utilized by NGS mapping/assembly applications.
At this time these file formats typically end with the following extensions:

    .ab1
    .sff
    .fastq
    .fasta
    .fna

ReadsBySample
=============

This directory contains only directories that are named after each of the samplenames that have been sequenced. The concept of this folder is to make it very easy to look up all data related to a given samplename.
