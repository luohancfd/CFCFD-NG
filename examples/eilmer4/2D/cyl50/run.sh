#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=cyl50
e4shared --run --job=cyl50 --verbosity=1 --max-cpus=4
e3post.py --job=cyl50 --tindx=all --vtk-xml
