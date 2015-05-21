#!/bin/bash
# run.sh
prep-gas ideal-air.inp ideal-air-gas-model.lua
e4shared --prep --job=rmi
e4shared --run --job=rmi --verbosity=1 --max-cpus=4 > LOGFILE_RUN
e3post.py --job=rmi --tindx=all --vtk-xml
