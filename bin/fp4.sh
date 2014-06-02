#!/bin/bash

#set -x

# wrapper for FP4 FP in Open Babel

if [ "$#" != "2" ] ; then
    echo "usage: "$0" query_molecule database"
    exit 1
fi

q=$1
db=$2

# score molecules
obabel $q $db -ofpt -xfFP4 2> /dev/null | grep Tanimoto | \
awk '/^>active/{print $6" 1"}/^>InChI/{print $6" 0"}' > $q".fp4.score-label"

# AUC and ER @ 1%
auc_tool -i $q".fp4.score-label" -p 0.01
