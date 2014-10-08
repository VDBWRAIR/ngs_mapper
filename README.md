The majority of the documentation is locate at:
https://vdbpm.org/projects/miseqpipeline/wiki

# Installation

1. Clone

  Assumes you already have git installed. If not you will need to get it installed by your system administrator.

  ```
  git clone https://github.com/VDBWRAIR/miseqpipeline.git
  cd miseqpipeline
  ```

2. Install System Packages

  This is the only part of the installation process that you should need to become the super user

  - Red Hat/CentOS(Requires the root password)
  
    ```
    su -c 'python setup.py install_system_packages'
    ```
  
  - Ubuntu
  
    ```
    sudo 'python setup.py install_system_packages'
    ```

3. Python

  The miseqpipeline requires python 2.7.3+ but < 3.0

  ```
  python setup.py install_python
  ```

  - Quick verify that Python is installed

    The following should return python 2.7.x(where x is somewhere from 3 to 9)

    ```
    $HOME/bin/python --version
    ```

4. Setup virtualenv
  
  
  1. Where do you want the pipeline to install? Don't forget this path, you will need it every time you want to activate the pipeline

    ```
    venvpath=$HOME/.miseqpipeline
    ```

  2. Install the virtualenv to the path you specified

    ```
    wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz#md5=f61cdd983d2c4e6aeabb70b1060d6f49 -O- | tar xzf -
    $HOME/bin/python virtualenv-1.11.6/virtualenv.py --prompt="(miseqpipeline) " $venvpath 
    ```

  3. Activate the virtualenv. You need to do this any time you want to start using the pipeline

    ```
    . $HOME/.miseqpipeline/bin/activate
    ```

5. Install the pipeline into virtualenv

  ```
  python setup.py install
  ```

  It should be safe to run this more than once in case some dependencies do not fully install.

6. Verify install

  You can pseudo test the installation of the pipeline by running the functional tests

  ```
  nosetests miseqpipeline/tests/test_functional.py
  ```

# Running the pipeline

  Before you use the pipeline you always need to ensure that you have the virtualenv activated that you installed into. Activating a virtualenv more than once is fine as it is smart enough to know if you already have done it. So if you are unsure, just do it anyways.
  
  If you copy pasted the installation instructions verbatim, then you can activate as follows:
  
  ```
  . $HOME/.miseqpipeline/bin/activate
  ```
  
  If you changed the venvpath in the installation then you will need to use that path instead of the above.

## Running a single sample

  You can run a single sample by using the runsample.py command. There are 2 examples that you can use. Sample 780 and Sample 947 which are both located in the
  miseqpipeline/tests/fixtures/functional directory.
  Inside that directory you will see directory for each sample which contains its reads as well as the reference for each of them and a config file for each sample. You can ignore the config file
  as it is used by the tests to determine if the sample ran correctly or not

  You can use runsample.py on them as follows:

  ```
  mkdir -p tdir && cd tdir
  runsample.py -od 780 ../miseqpipeline/tests/fixtures/functional/780{,.ref.fasta} 780
  runsample.py -od 947 ../miseqpipeline/tests/fixtures/functional/947{,.ref.fasta} 947
  ```

  This will create a temporary directory called tdir and cd into it then run both sample 780 as well as 947
  and put their results inside of their own directory named after themselves.

  From there you can explore them on your own
