#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=bump
e4shared --run --job=bump --max-cpus=4 --verbosity=1
e3post.py --job=bump --tindx=all --vtk-xml
