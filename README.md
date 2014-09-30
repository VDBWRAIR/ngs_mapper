The majority of the documentation is locate at:
https://vdbpm.org/miseqpipeline/Wiki

# Installation

1. Install System Packages

  This is the only part of the installation process that you should need to become the super user

  - Red Hat/CentOS
  
    ```
    su -c 'yum groupinstall "Development tools"; yum install -y wget ncurses{,-devel} zlib{,-devel} freetype{,-devel} readline{,-devel} openssl{,-devel} libpng{,-devel} ImageMagick java-1.7.0-openjdk git'
    ```
  
  - Ubuntu
  
    ```
    sudo apt-get install -y build-essential libncurses5{,-dev} zlib1g{,-dev} libpango1.0-{0,dev} libreadline6{,-dev} openssl libssl-dev unzip imagemagick libpng12-dev default-jre git
    ```

2. Clone

  ```
  git clone https://github.com/VDBWRAIR/miseqpipeline.git
  cd miseqpipeline
  ```

3. Python

  The miseqpipeline requires python 2.7.3+ but < 3.0

  1. Install Python 2.7.3+ into your home directory

    Ubuntu starting with Precise(12.04) has had Python 2.7.3+ so you technically 
    should not have to do this, but best to do it just to be safe.

    On RedHat/CentOS 6 and earlier this is mandatory since only python 2.6.6 comes
    shipped.

    ```
    mkdir -p ~/src && pushd ~/src
    prefix=$HOME
    version=2.7.8
    wget --no-check-certificate https://www.python.org/ftp/python/${version}/Python-${version}.tgz -O- | tar xzf -
    cd Python-${version}
    ./configure --prefix $prefix
    make
    make install
    popd
    ```

  2. Quick verify that Python is installed

    The following should return python 2.7.x(where x is somewhere from 3 to 9)

    ```
    $HOME/bin/python --version
    ```

4. Setup virtualenv

  1. Where do you want the pipeline to install? Don't forget this path, you will need it every time you want to activate the pipeline

  ```
  venvpath=/home/myusername-changeme/.miseqpipeline
  ```

  2. Install the virtualenv to the path you specified

  ```
  wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz#md5=f61cdd983d2c4e6aeabb70b1060d6f49 -O- | tar xzf -
  $HOME/bin/python virtualenv-1.11.6/virtualenv.py --prompt="(miseqpipeline) " $venvpath 
  ```

  3. Activate the virtualenv. You need to do this any time you want to start using the pipeline

  ```
  . /home/myuserename-changeme/bin/activate
  ```

5. Ensure pre-requisites

  Read the REQUIREMENTS.md file

6. Install into virtualenv

  ```
  python setup.py install
  ```

  It should be safe to run this more than once in case some dependencies do not fully install.

7. Verify install

  You can pseudo test the installation of the pipeline by running the functional tests

  ```
  nosetests miseqpipeline/tests/test_functional.py
  ```

# Running a single sample

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
