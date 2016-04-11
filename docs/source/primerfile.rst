Primer File
===========

You can set the primer file and the the primer options in the config.yaml file. You'll want to change the fields under trim_reads. To learn more about the options, see the help fields below, and see the documentation at http://www.usadellab.org/cms/?page=trimmomatic

.. code-block:: yaml

    trim_reads:
        primerfile:
            default:
            help:  'Primer File'
        primerseed:
            default: 2
            help:   'seed mismatches'
        palindromeclip:
            default: 30
            help: 'palindrome clip threshold'
        simpleclip:
            default: 20
            help: 'simple clip threshold' 

The available options are:

.. code-block:: yaml 

    primerfile: path to the primer file
    primerseed: seed mismatches
    palindromeclip: palindrome clip threshold
    simpleclip: simple clip threshold
