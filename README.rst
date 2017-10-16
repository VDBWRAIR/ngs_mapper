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

.. image:: https://img.shields.io/docker/build/tyghe/ngs_mapper.svg?style=plastic
    :target: https://hub.docker.com/r/tyghe/ngs_mapper

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

Running via Docker
^^^^^^^^^^^^^^^^^^

The following requires that you have docker installed on your computer.
You can get docker from visiting https://www.docker.com/ and clicking on the
Operating System from the 'Get Docker' dropdown.

Example on how to run runsample with Docker
+++++++++++++++++++++++++++++++++++++++++++

The following is simply an example of how to run the test data set that resides
in the ngs_mapper git repository. This is to simply get you familiar with
how to run the pipeline using docker.

To run the following example you will need to clone the ngs_mapper git repository
to your local computer

.. code-block:: bash

     git clone https://github.com/VDBWRAIR/ngs_mapper.git

Once you have cloned the repo you can then run through the following.

In the below command we are instructing docker to mount the 
`ngs_mapper/tests/fixtures/functional` directory from your local computer to 
`/NGSDATA` inside of the docker container when it runs.

It is also mounting the current directory you are running the command from as
`/output` inside the container.

You can then see that it is running the runsample command the same as you would
normally run it and pointing to the /data/947 directory to read samples data
from(as well as the reference file) and then sending all the output to the
/output/947 directory.

After the pipeline is finished running you would find the 947 directory under
the output directory in your local computer.

The `vdbwrair/ngs_mapper:latest` instructs docker to use the latest version
of ngs_mapper from dockerhub. 

If you want to use a different version you can replace `latest` with any tag
listed under https://hub.docker.com/r/vdbwrair/ngs_mapper/tags/

For example: `vdbwrair/ngs_mapper:v1.5.4` would use the `v1.5.4` version of
ngs_mapper.

.. code-block:: bash

    mkdir -p output
    docker run -it -v $PWD/output:/output -v $PWD/ngs_mapper/tests/fixtures/functional:/NGSDATA vdbwrair/ngs_mapper:latest runsample /NGSDATA/947 /NGSDATA/947.ref.fasta -od /output/947 947

Once you have completed this you should have a better understanding of how
to use docker to run ngs_mapper. To run your own data you can replace

.. code-block:: bash

    $PWD/ngs_mapper/tests/fixtures/functional

with the path to your data and then change the /NGSDATA/947 to reference
where your samples would reside.

For example, if your samples were in /some/path/sampledata/sample1, then 
you would use something like the following:

.. code-block:: bash

    mkdir -p output
    docker run -it -v $PWD/output:/output -v /some/path/sampledata:/NGSDATA vdbwrair/ngs_mapper:latest runsample /NGSDATA/sample1 /NGSDATA/ref.fasta -od /output/947 947

If your reference file exists somewhere outside of /some/path/sampledata you
can use another -v option for docker to make it available within the docker
container when it runs

.. code-block:: bash

    docker run -it -v /path/to/references:/references -v $PWD/output:/output -v /some/path/sampledata:/NGSDATA vdbwrair/ngs_mapper:latest runsample /NGSDATA/sample1 /references/ref.fasta -od /output/947 947
