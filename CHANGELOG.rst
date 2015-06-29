Changelog
---------

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
