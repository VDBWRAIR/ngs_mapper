Upgrade
=======

If you activate the pipeline via ngs_mapper/setup
----------------------------------------------------

  - Completely remove the existing ngs_mapper directory. 
  - Then follow :doc:`install`

If you installed using setup.py
-------------------------------------------

    #. First fetch any possible updates

        .. code-block:: bash
        
            cd ~/ngs_mapper; git fetch

    #. Then check if you need to update

        .. code-block:: bash

            git status | grep -q 'Your branch is behind' && echo 'You need to update' || echo 'You are up-to-date'
    
        If it returns You are up-to-date you are done

    #. Update(pull new code)

        .. code-block:: bash

            git pull

    #. Go into your ngs_mapper directory and rerun the setup script

        .. code-block:: bash

          python setup.py install

