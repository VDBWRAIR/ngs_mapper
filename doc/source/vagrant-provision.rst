====================
vagrant-provision.py
====================

.. py:module:: .vagrant-provision

Just a wrapper around the manual installation process.

Comes with two options:

* vagrant-provision.py --install-system-packages
* vagrant-provision.py --install-pipeline

Install System Packages
-----------------------

Just runs the following in a command shell

    .. code-block:: bash

        python setup.py install_system_packages

Install Pipeline
----------------

Essentially runs all the steps after the system package installation for you

* git clone to ~/miseqpipeline
* python setup.py install_python
* setup virtualenv into ~/.miseqpipeline
* activate ~/.miseqpipeline
* python setup.py install
