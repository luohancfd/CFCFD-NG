#!/bin/bash
# This test case gives an example EQ-spectra.py
# usage for shock heated air.
# Here we calculate the post-shock equilibrium 
# radiation intensity spectrum from 100 to 1200nm
# for shot 46-45 (see EAST.inp for the experiment 
# conditions).
# NOTE: The CEA program and scipy need to be installed
#       for this example to run.

# 1. Setup the radiation model LUA input file
script_rad2.py -i air-radiators.py -L rad-model.lua

# 2. Run EQ-spectra.py
EQ-spectra.py --input-file=EAST.py

# The files coefficient_spectra.txt and intensity_spectra.txt
# should now be in the working directory.
