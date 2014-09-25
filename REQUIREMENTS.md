# Requirements outside of python

- Python >=2.7.3
- samtools ==0.1.19[commit 96b5f2294ac0054230e88913c4983d548069ea4e]
- bwa ==0.7.6a

# Installation

  You may or may not be able to install these via your package manager(yum/apt)

  Typically you can download the tar.gz files for them and run the following commands
  
  ```
  tar xzf filename.tgz
  cd filename
  ./configure
  make
  make install
  ```

  Make sure you first create a src directory where we can download dependencies to and install them from

  ```
  mkdir src
  ```

## Python

  The miseqpipeline requires python 2.7.3+ but < 3.0

  - Red Hat

  1. Install System Packages

    ```
    su -c 'yum groupinstall "Development tools"'
    su -c "yum install -y yum install -y ncurses{,-devel} zlib{,-devel} freetype{,-devel} readline{,-devel} openssl{,-devel}"
    ```

  2. Install Python 2.7.3+ into your home directory

    ```
    prefix=$HOME
    version=2.7.8
    wget --no-check-certificate https://www.python.org/ftp/python/${version}/Python-${version}.tgz -O- | tar xzf -
    cd Python-${version}
    ./configure --prefix $prefix
    make
    ```
    Then copy the binary to somewhere in your PATH

  3. Quick verify that Python is installed

    The following should return python 2.7.x(where x is somewhere from 3 to 9)

    ```
    $HOME/bin/python --version
    ```

## Samtools

  ```
  pushd ~/src
  git clone https://github.com/samtools/samtools
  cd samtools
  git checkout 96b5f2294ac005423
  make
  eval cp -f $(pwd)/samtools $(dirs +1)/env/bin/samtools
  eval cp -f $(pwd)/bcftools/bcftools $(dirs +1)/env/bin/bcftools
  popd
  ```

## BWA

  ```
  pushd ~/src
  git clone https://github.com/lh3/bwa
  cd bwa
  git checkout 0.7.6a
  make
  eval cp -f $(pwd)/bwa $(dirs +1)/env/bin/bwa
  popd
  ```
