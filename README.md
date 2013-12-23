Undocumented scripts that are useful

* Requirements

** Normal
* Numpy
* Biopython
* pyBWA

** Testing
* Python Nose
* Python Mock

## Running run_bwa.py

* If you are providing the path to a directory of references you will want to
  compile them all into a single fasta file and index them before you run on a bunch of samples.
  run_bwa.py by default looks for a file named reference.fa in the current directory from which it was executed
  and will use it if the reference argument is not supplied. What this means is that if you make a project directory for each
  sample you map it will compile and index the references in the given directory for every sample which may slow the analysis
  down.
