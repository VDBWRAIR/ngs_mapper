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

Installation
============

1. Clone

    Assumes you already have git installed. If not you will need to get it installed by your system administrator.

    .. code-block:: bash

        git clone https://githubusername@github.com/VDBWRAIR/miseqpipeline.git
        cd miseqpipeline

.. _install-system-packages:

2. Install System Packages

    This is the only part of the installation process that you should need to become the super user

    - Red Hat/CentOS(Requires the root password)
  
        .. code-block:: bash

            su -c 'python setup.py install_system_packages && chmod 666 setuptools*'
  
    - Ubuntu
  
        .. code-block:: bash

            sudo 'python setup.py install_system_packages && chmod 666 setuptools*'

3. Configure the defaults

    You need to configure the miseqpipeline/config.yaml file.

    1. Copy the default config to config.yaml

        .. code-block:: bash

            cp miseqpipeline/config.yaml.default miseqpipeline/config.yaml

    2. Then edit the config file which is in `yaml <http://docs.ansible.com/YAMLSyntax.html>`_ format

        The most important thing is that you edit the NGSDATA value so that it contains the path to your NGSDATA directory

4. Python

    The miseqpipeline requires python 2.7.3+ but < 3.0

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

            venvpath=$HOME/.miseqpipeline

    2. Install the virtualenv to the path you specified

        .. code-block:: bash

            wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz#md5=f61cdd983d2c4e6aeabb70b1060d6f49 -O- | tar xzf -
            $HOME/bin/python virtualenv-1.11.6/virtualenv.py --prompt="(miseqpipeline) " $venvpath 

      3. Activate the virtualenv. You need to do this any time you want to start using the pipeline

            .. code-block:: bash

                . $HOME/.miseqpipeline/bin/activate

6. Install the pipeline into virtualenv

    .. code-block:: bash

        python setup.py install

    It should be safe to run this more than once in case some dependencies do not fully install.

7. Build and view complete documentation

    .. code-block:: bash

        cd doc
        make clean && make html
        firefox build/html/install.html
        cd ..

8. Verify install

    You can pseudo test the installation of the pipeline by running the functional tests

    .. code-block:: bash

        nosetests miseqpipeline/tests/test_functional.py
