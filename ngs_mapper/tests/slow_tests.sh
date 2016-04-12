#!/usr/bin/env bash

THIS=$(cd $(dirname ${BASH_SOURCE[0]}) && pwd)
nosetests -v $THIS/test_functional.py
