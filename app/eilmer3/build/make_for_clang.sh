#!/bin/bash
# script to set the compiler environment variables and then build.
# Usage:
# $ ./make_for_clang.sh 
# $ ./make_for_clang.sh install

export OMPI_CC=clang
export OMPI_CXX=clang
make TARGET=for_clang_openmpi $@
