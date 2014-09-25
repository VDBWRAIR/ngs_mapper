#!/usr/bin/env python

from miseqpipeline import run_bwa
import logging

logging.basicConfig(level=logging.DEBUG)

run_bwa.main()
