#!/bin/bash

#set -x

# output up top N first molecules of given MOL2 file
# (useful to cap the number of ligands per receptor to trim down
# a dataset)
# USAGE: mol2_head.sh NB_MOLECULES MOL2_FILE

nb_mols_asked=$1
input=$2

nb_mols_in_file=`egrep -c MOLECULE $input`

if [ $nb_mols_in_file -le $nb_mols_asked ]; then
    cat $input
else
    remove_after=$(( $nb_mols_asked + 1))
    after_last_line=`awk '{print NR":"$0}' $input | \
                     egrep MOLECULE | head -$remove_after | \
                     tail -1 | cut -d':' -f1`
    to_keep=$(( $after_last_line - 1 ))
    head -$to_keep $input
fi
