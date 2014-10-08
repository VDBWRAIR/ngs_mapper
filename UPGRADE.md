# If you installed using the install.sh or install_system.sh(aka you activate the pipeline via . miseqpipeline/setup)

  - Completely remove the existing miseqpipeline directory. 
  - Then follow "README.md":https://github.com/VDBWRAIR/miseqpipeline/blob/master/README.md file

# If you installed using (setup.py)

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
    
  5. Follow the install section in the (README.md)
