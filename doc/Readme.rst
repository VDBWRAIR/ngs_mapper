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
