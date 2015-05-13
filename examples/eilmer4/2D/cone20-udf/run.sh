#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cone20
e4shared --run --job=cone20 --verbosity=1 --max-cpus=2
e3post.py --job=cone20 --tindx=all --vtk-xml
