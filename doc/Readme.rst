Documentation is built using `sphinx <www.sphinx-doc.org>`_ using `restructuredText <http://sphinx-doc.org/rest.html>`_. There is a lot to learn, but once you get the hang of it it isn't too bad.

Building the documentation
--------------------------

.. code-block:: bash
   
    cd docs
    make clean && make html

Viewing the documentation
-------------------------

.. code-block:: bash

    firefox doc/build/html/index.html
