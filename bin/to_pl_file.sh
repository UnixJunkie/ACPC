#!/bin/bash

# convert a .mol2 to a .pl file
out=`echo $1 | sed 's/mol2$/pl/g'`

# this needs my patched version of Open Babel from
# https://github.com/UnixJunkie/openbabel/tree/logP_contrib_per_heavy_atom

~/usr/bin/obabel $1 -omol2 -O/dev/null --append logP | egrep '^BEGIN|^MOL |^ATOM |^END' > $out
