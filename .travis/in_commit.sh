#!/bin/bash

# Exit if not found
set -e

# Files that were modified in this commit
FILES_CHANGED=$(git diff --name-only HEAD~1)
echo "Files changed in current commit:"
echo ${FILES_CHANGED}

# Exit if files
for file in $@
do
   echo ${FILES_CHANGED} | grep -q $file || ( echo "$file not in commit"; exit 1 )
done
