===========
Development
===========

Contributing to the pipeline is fairly straight forward. The process is as follows:

#. Fork the VDBWRAIR `Miseqpipeline project <https://github.com/VDBWRAIR/miseqpipeline>`_ on GitHub
#. git clone your forked version to your local computer
#. Make changes to the code and ensure everything is tested
    * [[Make Tests]]
#. Once you have tested all your changes you should commit them and push them up to your github fork of the project
#. After you have pushed your changes to your fork on github you can create a pull request which essentially notifies the miseqpipeline maintainers that you have changes that you would like to apply and they can try them out.

Test Environment
================

The easiest way to ensure that the installer and everything works is to bring up a blank virtual machine to test inside of

The project is configured with a `Vagrant <https://www.vagrantup.com/>`_ file to make this easier

The Vagrant file that comes with the pipeline is configured to automatically :doc:`provision <vagrantfile>` either a CentOS 6.5 or Ubuntu 14.04 virtual machine.
You can bring either or both up with one of the following commands:

* CentOS 6.5

    .. code-block:: bash

        vagrant up centos65

* Ubuntu 14.04

    .. code-block:: bash

        vagrant up ubuntu1404
* Both

    .. code-block:: bash

        vagrant up

Documentation
=============

Documentation is built using `sphinx <www.sphinx-doc.org>`_ using `restructuredText <http://sphinx-doc.org/rest.html`_. There is a lot to learn, but once you get the hang of it it isn't too bad.

Building the documentation
--------------------------

.. code-block:: bash
   
    cd docs
    make clean && make html

Viewing the documentation
-------------------------

* Option 1:
    View on your computer

    .. code-block:: bash

        firefox docs/build/index.html

* Option 2:
    View on `github <https://github.com/necrolyte2/miseqpipeline/tree/v1.0/doc/source>`_ (not as nice)
