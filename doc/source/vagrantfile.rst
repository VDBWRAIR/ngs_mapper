===========
Vagrantfile
===========

This is the configuration file for `Vagrant <http://www.vagrantup.com>`_
It instructs vagrant on how to provision each of the virtual machines.

Essentially it brings up whichever vagrant box you tell it in the vagrant up command.
It ensures that any OS dependant things are handled prior to running :py:mod:`vagrant-provision`
When it is finished the Vagrant box you bring up should have the pipeline and all dependencies completely installed

CentOS 6.5
==========

You refer to this machine using centos65 in any of the vagrant commands:

.. code-block:: bash

    vagrant up centos65

Ubuntu 14.04
============

You refer to this machine using ubuntu1404 in any of the vagrant commands

.. code-block:: bash

    vagrant up ubuntu1404

How to use vagrant
==================

You will need to read the documentation at www.vagrantup.com
If you are just interested in how to bring up or destroy a vagrant box you can skip directly to the `Command Line Interface <http://docs.vagrantup.com/v2/cli/index.html>`_
