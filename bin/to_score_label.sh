#!/bin/bash

in=$1
out=$2

awk '/^active/{print $2" 1"}!/^active/{print $2" 0"}' $in > $out
