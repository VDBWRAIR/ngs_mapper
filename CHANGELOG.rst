Changelog
---------

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

Version 1.1.0
+++++++++++++

- Documentation updates
- Platforms now identified via identifiers inside read files instead of filenames
- IonTorrent sync added
- Various bug fixes
- base_caller.py can now utilize multiple processes to speed up analysis
- Documentation now installs with the pipeline
- run_bwa no longer makes temp directory but instead uses output path
