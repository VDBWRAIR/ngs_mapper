The majority of the documentation is locate at:
https://vdbpm.org/miseqpipeline/Wiki

# Installation

1. Clone

  ```
  git clone https://github.com/VDBWRAIR/miseqpipeline.git
  cd miseqpipeline
  ```

2. Python

  The miseqpipeline requires python 2.7.3+ but < 3.0

  - Red Hat

  1. Install System Packages

    ```
    su -c 'yum groupinstall "Development tools"'
    ```

    ```
    su -c "yum install -y yum install -y ncurses{,-devel} zlib{,-devel} freetype{,-devel} readline{,-devel} openssl{,-devel}"
    ```

  2. Install Python 2.7.3+ into your home directory

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

  3. Quick verify that Python is installed

    The following should return python 2.7.x(where x is somewhere from 3 to 9)

    ```
    $HOME/bin/python --version
    ```

3. Setup virtualenv

  ```
  wget --no-check-certificate https://pypi.python.org/packages/source/v/virtualenv/virtualenv-1.11.6.tar.gz#md5=f61cdd983d2c4e6aeabb70b1060d6f49 -O- | tar xzf -
  $HOME/bin/python virtualenv-1.11.6/virtualenv.py env 
  . env/bin/activate
  ```

4. Ensure pre-requisites

  Read the REQUIREMENTS.md file

5. Install into virtualenv

  ```
  python setup.py install
  ```

  It should be safe to run this more than once in case some dependencies do not fully install.
