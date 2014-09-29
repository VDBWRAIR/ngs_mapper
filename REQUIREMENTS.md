# Requirements outside of python

- Python >=2.7.3
- samtools ==0.1.19[commit 96b5f2294ac0054230e88913c4983d548069ea4e]
- bwa ==0.7.6a
- Trimmomatic ==0.32

# Installation

  You may or may not be able to install these via your package manager(yum/apt) but these instructions
  will help you manually install them.

  You need to ensure that you have a python virtualenv setup and activated prior to installing as these
  instructions assume that you do. The README.md has instructions on how to do that.

## System Packages

  - RedHat/CentOS

    ```
    su -c "yum install -y libpng{,-devel} ImageMagick java-1.7.0-openjdk"
    ```

  - Ubuntu

    ```
    sudo apt-get install -y imagemagick libpng12-dev default-jre
    ```

## Samtools

  ```
  mkdir -p ~/src && pushd ~/src
  git clone https://github.com/samtools/samtools
  cd samtools
  git checkout 96b5f2294ac005423
  make
  eval cp -f $(pwd)/samtools ${VIRTUAL_ENV}/bin/samtools
  eval cp -f $(pwd)/bcftools/bcftools ${VIRTUAL_ENV}/bin/bcftools
  popd
  ```
  - Verify Samtools

  ```
  samtools 2>&1 | grep -q 'Version: 0.1.18' && echo "samtools installed" || echo "Incorrect version or not installed"
  ```

## BWA

  ```
  mkdir -p ~/src && pushd ~/src
  git clone https://github.com/lh3/bwa
  cd bwa
  git checkout 0.7.6a
  make
  eval cp -f $(pwd)/bwa ${VIRTUAL_ENV}/bin/bwa
  popd
  ```

  - Verify BWA

  ```
  bwa 2>&1 | grep -q 'Version: 0.7.6a' && echo "BWA installed" || echo "Incorrect version or not installed"
  ```

## Trimmomatic

  ```
  mkdir -p ~/src && pushd ~/src
  wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
  unzip Trimmo*.zip
  cp -R Trimmo*/ ${VIRTUAL_ENV}/lib/
  popd
  ```

  - Verify Trimmomatic

  ```
  [[ $(java -jar /home/vmuser/miseqpipeline/env/lib/Trimmomatic-0.32/trimmomatic-0.32.jar 2>&1 | egrep 'SE|PE' | wc -l) == 2 ]] && echo "Trimmomatic installed" || echo "Trimmomatic not installed"
  ```
