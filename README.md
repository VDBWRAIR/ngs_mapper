# Installation

  [install](doc/source/install.rst)

# Upgrading
  
  [upgrade](doc/source/upgrade.rst)

# Changelog

  [changelog](CHANGELOG.rst)

# Running the pipeline

  Before you use the pipeline you always need to ensure that you have the virtualenv activated that you installed into. Activating a virtualenv more than once is fine as it is smart enough to know if you already have done it. So if you are unsure, just do it anyways.
  
  If you copy pasted the installation instructions verbatim, then you can activate as follows:
  
  ```
  . $HOME/.ngs_mapper/bin/activate
  ```
  
  If you changed the venvpath in the installation then you will need to use that path instead of the above.

## Running a single sample

  You can run a single sample by using the runsample.py command. There are 2 examples that you can use. Sample 780 and Sample 947 which are both located in the
  ngs_mapper/tests/fixtures/functional directory.
  Inside that directory you will see directory for each sample which contains its reads as well as the reference for each of them and a config file for each sample. You can ignore the config file
  as it is used by the tests to determine if the sample ran correctly or not

  You can use runsample.py on them as follows:

  ```
  mkdir -p tdir && cd tdir
  runsample.py -od 780 ../ngs_mapper/tests/fixtures/functional/780{,.ref.fasta} 780
  runsample.py -od 947 ../ngs_mapper/tests/fixtures/functional/947{,.ref.fasta} 947
  ```

  This will create a temporary directory called tdir and cd into it then run both sample 780 as well as 947
  and put their results inside of their own directory named after themselves.

  From there you can explore them on your own
