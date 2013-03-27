#!/bin/bash
# script to set the compiler environment variables and then build.
# Usage:
# $ ./make_for_gnu.sh 
# $ ./make_for_gnu.sh install

export OMPI_CC=gcc
export OMPI_CXX=g++
make TARGET=for_openmpi $@
