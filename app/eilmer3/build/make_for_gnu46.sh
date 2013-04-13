#!/bin/bash
# script to set the compiler environment variables and then build.
# Usage:
# $ ./make_for_gnu46.sh 
# $ ./make_for_gnu46.sh install

export OMPI_CC=gcc-4.6
export OMPI_CXX=g++-4.6
make TARGET=for_openmpi $@
