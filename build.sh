#!/bin/bash

#set -x

# build
obuild clean
obuild configure
obuild build
# install
ln -sf $PWD/dist/build/acpc/acpc                   ~/bin
ln -sf $PWD/dist/build/acpc_big/acpc_big           ~/bin
ln -sf $PWD/dist/build/acpc_mol2tool/acpc_mol2tool ~/bin
