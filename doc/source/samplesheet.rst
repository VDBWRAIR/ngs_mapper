Samplesheet
===========

You just need to make a file somewhere that is space or tab delimited with the first column being the sample name and the second column being the path to the reference to use to map that sample.
An easy way to get all the samplenames is by looking inside the MiSeq output directory where you will find a SampleSheet.csv file that you can open with excel/libreoffice

The resulting file should look very similar to this::

    my_precious_sample1 /path/to/references/reference1.fasta
    my_precious_sample2 reference2.fasta
    sample3 /path/to/contatenated_references.fasta

Note: This file does not contain any headers, but any line that starts with a # will be ignored so you could make the first line as follows if you want to::

    # Samplename Reference

Concatenate References
----------------------

Command to concatenate fasta files for more than one reference if you intend to try multiple references.

.. code-block:: bash

    cat file1.fasta file2.fasta file3.fasta > Flu__ConcatRefs.fasta

Note
----

Remember: If you are running you samples against a concatenated set of references, you will need to look at the output to determine which reference worked best for a given sample then re-run the pipeline on that sample or set of samples with the specific reference to ensure maximum use of data for that sample in assembly.

Generate Samplesheeet from MiSeq run
------------------------------------

Pulls out every samplename from a MiSeq SampleSheet.csv and creates a file that has samplename<TAB>REFPATH<NEWLINE>

    .. code-block:: bash

        awk -F',' '/^[0-9][0-9][0-9],/ {printf("%s\tREFPATH\n",$2);} /path/to/miseqrun/SampleSheet.csv > samplesheet.tsv

Command to generate a samplesheet from a 454 Run File
-----------------------------------------------------

Pulls out the column with the samplename and the column with the reference

    .. code-block:: bash

        awk '!/^[!#]/{printf("%s\t%s\n",$2,$6)}' runfile.txt > samplesheet.tsv
