# If you activate the pipeline via miseqpipeline/setup

  - Completely remove the existing miseqpipeline directory. 
  - Then follow [INSTALL.md](INSTALL.md)

# If you installed using [setup.py](setup.py)

  1. First fetch any possible updates

    ```cd ~/miseqpipeline; git fetch```

  2. Then check if you need to update

    ```
    git status | grep -q 'Your branch is behind' && echo 'You need to update' || echo 'You are up-to-date'
    ```
    
    - If it returns You are up-to-date you are done

  3. Update(pull new code)

    ```
    git pull
    ```
    
  4. Ensure your virtualenv is activated(same as before you run the pipeline)

    ```
    . ~/.miseqpipeline/bin/activate
    ```
    
  5. Run setup.py install

    ```
    python setup.py install
    ```
