#!/bin/bash
# This simple test case demonstrates how to use the 
# EQ_spectra.py program.
# Here we calculate the post-shock equilibrium 
# radiation intensity spectrum from 100 to 2000nm
# for a hypothetical Argon shock tube condition.
# NOTE: The CEA program and scipy need to be installed
#       for this example to run.

# 1. Setup the radiation model LUA input file
radmodel.py -i argon-radiators.py -L rad-model.lua

# 2. Run EQ_spectra.py
EQ_spectra.py --input-file=argon.py

# The files coefficient_spectra.txt and intensity_spectra.txt
# should now be in the working directory.
