#!/bin/bash

# list molecule names in MOL2 files

egrep --no-filename -A1 MOLECULE "$@" | egrep -v MOLECULE | grep -v '\-\-'
