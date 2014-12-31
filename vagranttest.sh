#!/bin/bash

# Simple script to test a given vagrant box

if [ -z "$1" ]
then
    echo "Must supply a boxname listed within the vagrant status command output"
    exit 1
fi

if vagrant status $1 | grep -q 'running'
then
    vagrant destroy $1
fi

# Bring up box. Provision and run tests
vagrant up $1
vagrant ssh $1 -c ". ~/.ngs_mapper/bin/activate; cd ngs_mapper; nosetests ngs_mapper"
