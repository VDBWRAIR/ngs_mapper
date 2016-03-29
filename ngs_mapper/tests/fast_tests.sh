#!/usr/bin/env bash

THIS=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
nosetests -v --nologcapture -a '!slow' -I '.*test_functional.py$' $THIS
