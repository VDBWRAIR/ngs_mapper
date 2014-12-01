============
Syncing Data
============

Sync User
=========

To ensure data is protected from accidental deletion/modification it is suggested that you :doc:`setup a sync user <syncuser>` on your system and run the sync scripts as that user.
This way if you accidentally run a command such as 

    .. code-block:: bash

        rm -rf /path/to/NGSData

You will only get permission denied instead of accidentally deleting your data

Running a sync command as your sync user
========================================

You will have to login to the sync user before running the commands in order to use that account to sync the data
This could mean that you use the su, sudo or just ssh as that user. This will depend on your system setup

The following scripts exist to sync data from an instrument into the data structure:

* :py:mod:`roche_sync <miseqpipeline.roche_sync>`
* :py:mod:`miseq_sync <miseqpipeline.miseq_sync>`
* :py:mod:`sanger_sync <miseqpipeline.sanger_sync>`
* iontorrent_sync.py coming soon
