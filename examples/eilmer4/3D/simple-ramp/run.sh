#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=ramp
e4shared --run --job=ramp --verbosity=1
e3post.py --job=ramp --tindx=all --vtk-xml
