#!/bin/bash

#set -x

# build
obuild configure
obuild build

# soft link .annot files so that Emacs' tuareg-mode can find them
mkdir -p _build
for f in `find dist | egrep '\.annot$'` ; do
    ln -sf ../$f _build/
done

# install
ln -sf $PWD/dist/build/acpc_pltool/acpc_pltool ~/bin
ln -sf $PWD/dist/build/acpc_pqrtool/acpc_pqrtool ~/bin
ln -sf $PWD/dist/build/acpc_mol2tool/acpc_mol2tool ~/bin
ln -sf $PWD/dist/build/acpc_codec/acpc_codec ~/bin
ln -sf $PWD/dist/build/acpc_mol2reader/acpc_mol2reader ~/bin
ln -sf $PWD/dist/build/acpc_auctool/acpc_auctool ~/bin
ln -sf $PWD/dist/build/acpc_ertool/acpc_ertool ~/bin
ln -sf $PWD/dist/build/acpc_consrank/acpc_consrank ~/bin
ln -sf $PWD/dist/build/acpc/acpc ~/bin
ln -sf $PWD/dist/build/acpc_scorer/acpc_scorer ~/bin
ln -sf $PWD/dist/build/acpc_par/acpc_par ~/bin

## helper scripts for research experiments, OB FP4 and MACCS wrappers, etc.
# ln -sf $PWD/bin/to_pl_file.sh     ~/bin
# ln -sf $PWD/bin/to_score_label.sh ~/bin
# ln -sf $PWD/bin/maccs.sh          ~/bin
# ln -sf $PWD/bin/fp4.sh            ~/bin
