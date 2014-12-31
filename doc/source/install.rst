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
    * At least 1GB per CPU core for small genomes(Dengue, Flu,...)

Python Packages
---------------

All python packages can be defined in a pip requirements.txt file
The pipeline comes with all of the necessary python packages already defined inside of `requirements.txt`_.

.. _requirements.txt: ../../../requirements.txt

System Packages
---------------

The pipeline requires some system level packages(software installed via your Linux distribution's package manager)
The installer looks for the `system_packages.lst <../../../system_packages.lst>`_ and installs the correct packages using that file.
This file is a simple json formatted file that defines packages for each package manager

Roche Utilities
^^^^^^^^^^^^^^^

If you intend on using the :py:mod:`roche_sync <ngs_mapper.roche_sync>` you will need to ensure that the ``sfffile`` command is in your PATH. That is, if you execute ``$> sfffile`` it returns the help message for the command.

This command should automatically be installed and put in your path if you install the Data Analysis CD #3 that was given to you with your Roche instrument.

MidParse.conf
^^^^^^^^^^^^^

If you inted on using the :py:mod:`roche_sync <ngs_mapper.roche_sync>` you may need to edit the included ngs_mapper/MidParse.conf file before installing. This file is formatted to be used by the Roche utilities and more information about how it is used can be found in the Roche documentation.

Installation
============

1. Clone

    Assumes you already have git installed. If not you will need to get it installed by your system administrator.

    .. code-block:: bash

        git clone https://githubusername@github.com/VDBWRAIR/ngs_mapper.git
        cd ngs_mapper

.. _install-system-packages:

2. Install System Packages

    This is the only part of the installation process that you should need to become the super user

    - Red Hat/CentOS(Requires the root password)
  
        .. code-block:: bash

            su -c 'python setup.py install_system_packages && chmod 666 setuptools*'
  
    - Ubuntu
  
        .. code-block:: bash

            sudo python setup.py install_system_packages && sudo chmod 666 setuptools*

3. Configure the defaults

    You need to configure the ngs_mapper/:doc:`config` file.

    1. Copy the default config to config.yaml

        .. code-block:: bash

            cp ngs_mapper/config.yaml.default ngs_mapper/config.yaml

    2. Then edit the ngs_mapper/config.yaml file which is in `yaml <http://docs.ansible.com/YAMLSyntax.html>`_ format

        The most important thing is that you edit the NGSDATA value so that it contains the path to your NGSDATA directory.

        **The path you use for NGSDATA must already exist**

        .. code-block:: bash

            mkdir -p /path/to/NGSDATA

4. Python

    The ngs_mapper requires python 2.7.3+ but < 3.0

    - Ensure python is installed

        .. code-block:: bash

            python setup.py install_python

    - Quick verify that Python is installed

        The following should return python 2.7.x(where x is somewhere from 3 to 9)

        .. code-block:: bash

            $HOME/bin/python --version

5. Setup virtualenv
  
  
    1. Where do you want the pipeline to install? Don't forget this path, you will need it every time you want to activate the pipeline

        .. code-block:: bash

            venvpath=$HOME/.ngs_mapper

    2. Install the virtualenv to the path you specified

        .. code-block:: bash

            wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz#md5=f61cdd983d2c4e6aeabb70b1060d6f49 -O- | tar xzf -
            $HOME/bin/python virtualenv-1.11.6/virtualenv.py --prompt="(ngs_mapper) " $venvpath 

      3. Activate the virtualenv. You need to do this any time you want to start using the pipeline

            .. code-block:: bash

                . $HOME/.ngs_mapper/bin/activate

6. Install the pipeline into virtualenv

    .. code-block:: bash

        python setup.py install

    It should be safe to run this more than once in case some dependencies do not fully install.


Build and view complete documentation
-------------------------------------

.. code-block:: bash

    cd doc
    make clean && make html
    firefox build/html/install.html#build-and-view-complete-documentation
    cd ..

Verify install
--------------

You can pseudo test the installation of the pipeline by running the functional tests

.. code-block:: bash

    nosetests ngs_mapper/tests/test_functional.py
