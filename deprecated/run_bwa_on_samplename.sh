#!/usr/bin/env bash

# Activate the precious
. /home/EIDRUdata/programs/vdbapps/bin/activate

# The samplename to work with
# If you change this then it should work on any samplename
samplename="$1"
if [ -z "$samplename" ]
then
    echo "I require a samplename to send my minion BWA to feast on"
    exit 1
fi

# What reference should be used?
reference="$2"
if [ -z "$reference" ]
then
    # No reference was specified so now we nail them with the huge list
    # Prompt the user with this every time
    PS3="Select a reference[quit to exit the menu]: "
    # Bash select statements are like pulling teeth until you figure them out
    select ref in /home/EIDRUdata/Analysis/References/*.fasta
    do
        # Make sure the user didn't type something stupid...they do that sometimes
        case $REPLY in
            "quit")
                if [ ! -z "$reference" ]
                then
                    break
                else
                    echo "You need to select a reference before quitting silly."
                fi
                ;;
            *)
                reference=$ref
                echo "You selected $ref"
                ;;
        esac
    done
fi
# Arr, where be the sample files?
readfilepath="/home/EIDRUdata/NGSData/ReadsBySample"

# Some bash expansion array garbage
readfiles=( ${readfilepath}/${samplename}/* )

# Ohh noes, what if there are other read files other than just 2 of them in the folder
# Superman please help me
if [ ${#readfiles[@]} -ne 2 ]
then
    # No silly, we will just bomb out because this is programming for biology
    # Cryptic error message as per the standard bio-programming
    echo "Not exactly 2 miseq read files found"
    exit 1
fi

# How many cpus can we thrash
cpus=$(cat /proc/cpuinfo | grep processor | wc -l)
# Don't use them all silly, let one be available for that game of minesweeper while we wait for BWA to finish
cpus=$(($cpus - 1))

# Yay, I think we can run
map_bwa.py -t $cpus $reference ${readfilepath}/${samplename}/${samplename}_*_R{1,2}*.fastq --output ${samplename}.sai
# Convert that icky tab file to a useful indexed compressed bam file
sai_to_bam ${samplename}.sai && rm ${samplename}.sai
