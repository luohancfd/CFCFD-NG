#!/bin/bash
# This test case gives an example EQ_spectra.py
# usage for shock heated air.
# Here we calculate the post-shock equilibrium 
# radiation intensity spectrum from 100 to 1200nm
# for shot 46-45 (see EAST.inp for the experiment 
# conditions).
#
# Two spectral models are considered:
# 1) Legacy/baseline spectral modelling (circa Potter PhD)
# 2) Advanced spectral modelling (circa 2013/2014)
# See sections 8.1 and 8.2 of the user-guide for a 
# description of each of these modellings.
#
# NOTE: The CEA program and scipy need to be installed
#       for this example to run.

# 1. Legacy/baseline spectral modelling
# 1a. Setup the radiation model LUA input file
radmodel.py -i air-radiators-legacy.py -L rad-model.lua

# 1b. Run EQ_spectra.py
EQ_spectra.py --input-file=EAST.py

# 1c. move output files 
mv coefficient_spectra.txt legacy_coefficient_spectra.txt
mv intensity_spectra.txt legacy_intensity_spectra.txt

# 2. New/advanced spectral modelling
# 2a. Setup the radiation model LUA input file
radmodel.py -i air-radiators-NIST-TB.py -L rad-model.lua

# 2b. Run EQ_spectra.py
EQ_spectra.py --input-file=EAST.py

# 2c. move output files 
mv coefficient_spectra.txt new_coefficient_spectra.txt
mv intensity_spectra.txt new_intensity_spectra.txt

# 3. Plot the result
gnuplot plot.gnu
