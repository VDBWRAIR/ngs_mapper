Changelog
---------

Version X
+++++++++

- Continuous Delivery support added for travis
- nfilter will now simply symlink if no options are supplied essentially skipping
  itself
- nfilter utilizes threads from config file
- config file now has THREADS default
- fix for bug where some miseq reads were not identified correctly in tagreads
- convert functions now support output directory
- bug fix for nfilter symlinking
- fix for qsub job output from runsample

Version 1.4.2
+++++++++++++

- Pipeline now handles gzip(.gz) input files
- Pipeline now handles ab1 input files
- Added Zenodo badge

Version 1.4.1
+++++++++++++

- IGV is installed with pipeline
- samtools version reverted back to same version as pre-1.4.0

Version 1.4.0
+++++++++++++

- Installation now utilizes miniconda to handle system dependencies such as
  bwa, samtools, trimmomatic, imagemagick. This is a substantial difference and will
  require a complete reinstall of the pipeline to upgrade.
  Miniconda installation removes a lot of code that needed to be maintained and
  streamlines the installation and makes it much faster.
- Added install.sh that makes installing/upgrading much easier.
  The tests also use this so the installation is tested much better now.
- Pipeline utilizes requirements-conda.txt to determine python+system software
  dependencies. This allows specifying versions and removes need for a
  system administrator to install.
- runsample now supports --primer-file option and other primer trimming options
  which will utilize trimmomatic's ILLUMINACLIP option
- runsamplesheet.sh supports an optional additional column in a given samplesheet
  that represents the primer fasta file to use to find sequences to trim out.
- Pipeline now looks for amount of threads instead of cpu cores. This will mean that
  on systems with hyperthreading that 2x more samples will run in parallel than before.
- Fixed bug where some parts of pipeline were not logging at all
- Fixed bug where graphs.sh could fail, yet pipeline would continue as if nothing
  was wrong
- Updated functional tests to include primer test
- Updated functional tests to output more information

Version 1.3.0
+++++++++++++

- Added ngs_filter stage/script that can filter based on index fastq files as well
  as reads that contain an N. This stage is off by default.
- Fixed a bug where some scripts were not logging properly

Version 1.2.4
+++++++++++++

- Fixes documentation issue with umask for sync user

Version 1.2.3
+++++++++++++

- Added travis-ci support to automatically run tests when code is pushed to github
- Projects now default to running inside of a temporary directory inside of the
  specified output directory(-od)
- runsample now sets TMPDIR to tmpdir inside of output directory so that all
  analysis is run within that directory  

Version 1.2.2
+++++++++++++

- runsample accepts --qsub_m and --qsub_l commands which will direct it to
  return a PBS qsub job that can be piped into qsub
- Added Python 2.6 support

Version 1.2.1
+++++++++++++

- Removed all occurances of bqd.mpileup and replaced with samtools.mpileup
- Changed bqd.parse_pileup such that it utilizes samtools.MPileupColumn to
  generate the dictionary items
- Remove legacy BamCoverage code that is not used anywhere
- Added support to select reads by specific platforms in runsample.py
- Fixed bug where MiSeq Index reads were being included in the mapping
- Renamed unpaired read file name that is produced by trim_reads from
  a generic Roche454 read name to simply unpaired_trimmed.fastq

Version 1.2.0
+++++++++++++

- Added reflen to qualdepth.json files since length only told you the length
  of the assembly and not the reference.
- Fixed issue where coverage graphic was not drawing gap lines at the end of
  references because there was no data.
- sample_coverage colors were hard to distinquish so they were changed
- Bug with sample_coverage where certain combinations of # of references
  and # of samples would generate a graphic where sub-plots for each reference
  were overlapping
- Fixed incorrect command in doc/README.rst for how to open documentation with Firefox
- Fixed issue with sample_coverage's usage statement and arguments description
- Fixed issue when no reads mapped and graphsample.py would raise an exception
- Fixed an issue when there were directories inside of the path specified that
  contains read files
- Replaced all .py scripts with same name but without .py. This is the correct
  way to have binary scripts for python. Aka, runsample.py is now just
  runsample

Version 1.1.0
+++++++++++++

- Documentation updates
- Platforms now identified via identifiers inside read files instead of filenames
- IonTorrent sync added
- Various bug fixes
- base_caller.py can now utilize multiple processes to speed up analysis
- Documentation now installs with the pipeline
- run_bwa no longer makes temp directory but instead uses output path
