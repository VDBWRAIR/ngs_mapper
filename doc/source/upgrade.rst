Upgrade
=======

If you activate the pipeline via miseqpipeline/setup
----------------------------------------------------

  - Completely remove the existing miseqpipeline directory. 
  - Then follow :doc:`install`

If you installed using setup.py
-------------------------------------------

    1. First fetch any possible updates

        .. code-block:: bash
        
            cd ~/miseqpipeline; git fetch

    2. Then check if you need to update

        .. code-block:: bash

            git status | grep -q 'Your branch is behind' && echo 'You need to update' || echo 'You are up-to-date'
    
        If it returns You are up-to-date you are done

    3. Update(pull new code)

        .. code-block:: bash

            git pull

    4. Reinstall

        Then follow :ref:`install <install-system-packages>`
