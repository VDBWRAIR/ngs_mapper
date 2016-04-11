=======
Install
=======

Requirements
============

Hardware
--------

* CPU
    * Quad Core 2.5GHz or better
        * More cores = faster run time when running multiple samples
        * Faster GHz = faster each sample runs
* RAM
    * This really depends on your data size
    
      If you are analyzing a 96 sample run then you should be fine with 1GB per CPU core
      
      If you are analyzing a 24 sample run then you will probably need about 4GB per CPU core since there will be more data

Roche Utilities
^^^^^^^^^^^^^^^

If you intend on using the :py:mod:`roche_sync <ngs_mapper.roche_sync>` you will 
need to ensure that the ``sfffile`` command is in your PATH. That is, if you 
execute ``$> sfffile`` it returns the help message for the command.

This command should automatically be installed and put in your path if you install 
the Data Analysis CD #3 that was given to you with your Roche instrument.

MidParse.conf
^^^^^^^^^^^^^

If you inted on using the :py:mod:`roche_sync <ngs_mapper.roche_sync>` you may need 
to edit the included ngs_mapper/MidParse.conf file before installing. This file is 
formatted to be used by the Roche utilities and more information about how it is 
used can be found in the Roche documentation.

Installation
============

1. Clone/Download the version you want

   #. Clone the code

      .. code-block:: bash

        git clone https://github.com/VDBWRAIR/ngs_mapper.git
        cd ngs_mapper
        
   #. Check which versions are available
   
      .. code-block:: bash
      
         git tag
   
   #. Checkout the version you want(current version |release|)
   
      .. code-block:: bash
      
         git checkout -b vX.Y.Z vX.Y.Z

2. Configure the defaults

    You need to configure the ngs_mapper/:doc:`config` file.

    1. Copy the default config to config.yaml

        .. code-block:: bash

            cp ngs_mapper/config.yaml.default ngs_mapper/config.yaml

    2. Then edit the ngs_mapper/config.yaml file which is in 
       `yaml <http://docs.ansible.com/YAMLSyntax.html>`_ format

        The most important thing is that you edit the NGSDATA value so that it 
        contains the path to your NGSDATA directory.

        **The path you use for NGSDATA must already exist**

        .. code-block:: bash

            mkdir -p /path/to/NGSDATA

3. Install

    The project now comes with a much more simplified installer which is based
    on miniconda.

    The following will install the project into the current directory that you
    are in.

    .. code-block:: bash

        ./install.sh miniconda

4. PATH Setup

    Once the project is installed you will need to setup your PATH environmental
    variable to include the 

    .. code-block:: bash

        export PATH=$PWD/miniconda/bin:$PATH

    You can put this into your .bashrc file inside your home directory so that any
    time you open a new terminal it automatically is run.
    
    If you don't setup your .bashrc you will have to run the export command from
    above every time you open a new terminal.

Verify install
--------------

You can pseudo test the installation of the pipeline by running the functional tests

.. code-block:: bash

    ngs_mapper/tests/slow_tests.sh

Documentation
-------------
The documentation is available to view online at http://ngs-mapper.readthedocs.org/en/latest/ By default this site always shows the latest documentation; you can select your verion by clicking `v:latest` in the bottom left menu, then selecting your version number from `Versions.`

If for any reason you need to use the documentation locally, you can build it:

.. code-block:: bash

    cd doc
    make clean && make html
    firefox build/html/install.html#build-and-view-complete-documentation
    cd ..
