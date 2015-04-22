#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=sharp
e4shared --run --job=sharp --verbosity=1
e3post.py --job=sharp --tindx=all --vtk-xml
