=========
Sync User
=========

Check to see if there is a sync user exists on the computer you are on
======================================================================

Assuming that the sync user was created with the name "ngssync" you can check for it as follows:

    .. code-block:: bash

        id ngssync

**Outputs either:**

#. User does not exist

    .. code-block:: bash

        id: ngssync: No such user

#. User exists

    .. code-block:: bash

        uid=500(ngssync) gid=500(ngssync) groups=500(ngssync)

    Which tells you that the user exists and the user id is 500, primary group id is 500, and a full list of all the groups that user is in. In this case ngssync is only in the one group wich is also the same as the primary group(ngssync)

Create the Sync User
====================

#. You need to be root to create a new user

    .. code-block:: bash

        su -
#. Setup some useful variables for later

    .. code-block:: bash

        sync_username="ngssync" # Change this if you want
        ngs_data_path="/path/to/your/NGSData" # Set this to your path

#. Then create the user account

    .. code-block:: bash

        useradd -s /bin/bash ${sync_username}
        passwd ${sync_username}

#. Setup the default umask for this account so files/directories are read only

    .. code-block:: bash

        echo "umask 0022" >> ~/.bash_profile

Ensure your sync user has read/write/execute to your data
=========================================================

Check access to files
---------------------

You need to make sure that your sync user has rwx access to your :doc:`NGS Data Structure <ngsdata>`

    .. code-block:: bash

        $> ls -l /path/to/NGSData
        dr-xr-xr-x.    8 ngssync ngssync   4096 Jul 14 11:56 RawData
        drwxr-xr-x.    7 ngssync ngssync   4096 Jul  9 08:07 ReadData
        drwxr-xr-x. 3667 ngssync root 286720 Sep 15 16:22 ReadsBySample

Listing format:

===============  ======  =========  =====  ====  ====  ==============
**permissions**  ignore  **owner**  group  size  date  directory/file
===============  ======  =========  =====  ====  ====  ==============

The important thing is that the ngssync is the username listed for the owner(3rd column) and that the permissions start with

    .. code-block:: bash

        drwx

So in the example above you can see that the RawData is owned by ngssync but does not have w(write) permissions which is wrong. You may think that ReadsBySample is setup incorrectly since the root group is set as the group(4th column), but that really doesn't matter since ngssync has rwx(read/write/execute) permissions.

If you see a directory that looks incorrect you can fix it with the following command:

    .. code-block:: bash

        su -c "chown -R ngssync /path/to/directory; chmod -R u=rwX,go=rX /path/to/directory"

**Requires root privileges**

Be patient as this can take a really long time if you have lots of data
