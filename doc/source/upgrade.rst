Upgrade
=======

If you activate the pipeline via miseqpipeline/setup
----------------------------------------------------

  - Completely remove the existing miseqpipeline directory. 
  - Then follow :doc:`install`

If you installed using setup.py
-------------------------------------------

    #. First fetch any possible updates

        .. code-block:: bash
        
            cd ~/miseqpipeline; git fetch

    #. Then check if you need to update

        .. code-block:: bash

            git status | grep -q 'Your branch is behind' && echo 'You need to update' || echo 'You are up-to-date'
    
        If it returns You are up-to-date you are done

    #. Update(pull new code)

        .. code-block:: bash

            git pull

    #. Reinstall

        Then follow :ref:`install <install-system-packages>`
