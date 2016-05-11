#!/bin/bash

# count molecules in MOL2 files

egrep -c MOLECULE "$@"
