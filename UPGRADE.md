# If you activate the pipeline via ngs_mapper/setup

  - Completely remove the existing ngs_mapper directory. 
  - Then follow [INSTALL.md](INSTALL.md)

# If you installed using [setup.py](setup.py)

  1. First fetch any possible updates

    ```cd ~/ngs_mapper; git fetch```

  2. Then check if you need to update

    ```
    git status | grep -q 'Your branch is behind' && echo 'You need to update' || echo 'You are up-to-date'
    ```
    
    - If it returns You are up-to-date you are done

  3. Update(pull new code)

    ```
    git pull
    ```

  4. Reinstall

    Then follow [INSTALL.md](INSTALL.md) starting with section #2
