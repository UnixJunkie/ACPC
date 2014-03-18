#!/bin/bash

#set -x

# configure, build and install via OASIS

oasis setup
ocaml setup.ml -configure --prefix ${HOME}/usr
ocaml setup.ml -build
ocaml setup.ml -install
