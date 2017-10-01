#!/usr/bin/env bash

#docker build -t ngs_mapper:1.5 .
docker run -it -v /Users/z0023fp/Projects/necrolyte2/ngs_mapper/NGSData:/NGSDATA ngs_mapper:1.5 ngs_mapper/tests/slow_tests.sh
