#!/bin/bash
# Author: Rowan J. Gollan
# Date: 08-Mar-2012
# Place: Bray Park, Queensland, Australia
#
# This script is used to parse the eilmer output
# for reporting norms and collecting it in a file.
#

# 1. Run eilmer to compute and print norm values
e3post.py --job=odw --tindx=all --gmodel-file=binary-gas.lua \
    --ref-function=odw-ref-function.py \
    --global-norm-list="rho,L1" > e3post.out


# 2. Run awk to collect the time and corresponding norms
#    Pipe to uniq to filter out repeated final tindx
awk -f parse-L1-norms.awk < e3post.out | uniq > rho-L1-norms.dat

