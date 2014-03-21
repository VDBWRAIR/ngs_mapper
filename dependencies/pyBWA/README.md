# pyBWA

Python libraries to interact with BWA mapper

The goal is to make it very easy to run bwa from any python script

## Install

```bash
git clone https://github.com/VDBWRAIR/pyBWA.git
pip install biopython
python setup.py install
```

## Simple Example

```python
import bwa
import sys

reference_path = sys.argv[1]
read_path = sys.argv[2]

# Ensure reference is indexed
bwa.index_ref( reference_path )

# Setup and run bwa mem
mem = bwa.BWAMem( reference_path, read_path )
retstat = mem.run( 'myoutput.sai' )

# Check return status
if retstat != 0:
    sys.stderr.write( "Error running bwa" )
```

## Multiple reads and references

```python
import bwa
import sys

reference_path = sys.argv[1]
read_path = sys.argv[2]

# If read_path or reference_path are a directory
# they can contain multiple files that will be concatted
# together into a single file.
# Reads can also be .sff files that will be converted to fastq
reads = bwa.compile_reads( read_path )
refs = bwa.compile_refs( reference_path )

# Ensure reference is indexed
bwa.index_ref( refs )

# Setup and run bwa mem
mem = bwa.BWAMem( refs, reads )
retstat = mem.run( 'myoutput.sai' )

# Check return status
if retstat != 0:
    sys.stderr.write( "Error running bwa" )
```

## Executables

pyBWA comes with some utility executables that wrap the functionality of BWA mapping
in a single easy to use executable.

* map_bwa.py wraps up bwa mem mapping in an easy to use single executable
  * It utilizes bwa.index_ref, bwa.compile_reads and bwa.compile_refs to easily 
    add SFF files, fastq files and fasta reference files to the mapping
* sai_to_bam converts the output sai sam file to an indexed/sorted bam file
