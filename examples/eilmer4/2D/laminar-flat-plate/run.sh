#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=lam_flat_plate_march
e4shared --run --job=lam_flat_plate_march --verbosity=1 --max-cpus=4
e3post.py --job=lam_flat_plate_march --tindx=all --vtk-xml
