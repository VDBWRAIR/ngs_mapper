.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.46716.svg
   :target: http://dx.doi.org/10.5281/zenodo.46716

.. image:: https://badge.waffle.io/VDBWRAIR/ngs_mapper.png?label=ready&title=Ready 
    :target: https://waffle.io/VDBWRAIR/ngs_mapper
    :alt: 'Stories in Ready'

.. image:: https://readthedocs.org/projects/ngs_mapper/badge/?version=latest
    :target: http://ngs_mapper.readthedocs.org/en/latest/
    :alt: Documentation Status

.. image:: https://travis-ci.org/VDBWRAIR/ngs_mapper.svg
    :target: https://travis-ci.org/VDBWRAIR/ngs_mapper

.. image:: https://coveralls.io/repos/VDBWRAIR/ngs_mapper/badge.svg
    :target: https://coveralls.io/r/VDBWRAIR/ngs_mapper

Installation
------------


`install <doc/source/install.rst>`_

Upgrading
---------
  
`upgrade <doc/source/upgrade.rst>`_

Changelog
---------

`changelog <CHANGELOG.rst>`_

Running the pipeline
--------------------

Before you use the pipeline you always need to ensure that the miniconda environment
is in your environment's path. You will want to look at the
`install <doc/source/install.rst>`_ for more details


Running a single sample
^^^^^^^^^^^^^^^^^^^^^^^

You can run a single sample by using the runsample.py command. There are 2 examples that you can use. Sample 780 and Sample 947 which are both located in the
ngs_mapper/tests/fixtures/functional directory.
Inside that directory you will see directory for each sample which contains its reads as well as the reference for each of them and a config file for each sample. You can ignore the config file
as it is used by the tests to determine if the sample ran correctly or not

You can use runsample.py on them as follows:

.. code-block:: bash

    mkdir -p tdir && cd tdir
    runsample -od 780 ../ngs_mapper/tests/fixtures/functional/780{,.ref.fasta} 780
    runsample -od 947 ../ngs_mapper/tests/fixtures/functional/947{,.ref.fasta} 947

This will create a temporary directory called tdir and cd into it then run both sample 780 as well as 947
and put their results inside of their own directory named after themselves.

From there you can explore them on your own
