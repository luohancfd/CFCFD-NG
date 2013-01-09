#!/bin/bash
# This simple test case demonstrates how to use the 
# EQ-spectra.py program.
# Here we calculate the post-shock equilibrium 
# radiation intensity spectrum from 100 to 2000nm
# for a hypothetical Argon shock tube condition.
# NOTE: The CEA program and scipy need to be installed
#       for this example to run.

# 1. Setup the radiation model LUA input file
script_rad2.py -i argon-radiators.py -L rad-model.lua

# 2. Run EQ-spectra.py
EQ-spectra.py --input-file=argon.py

# The files coefficient_spectra.txt and intensity_spectra.txt
# should now be in the working directory.
