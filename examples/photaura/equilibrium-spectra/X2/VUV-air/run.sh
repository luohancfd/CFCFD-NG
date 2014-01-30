#!/bin/bash
# This test case is from Umar Sheikh's thesis, where he analysed
# the steady region of the shock layer in front of a flat rectangular
# model (as close to 'equilibrium' as we can get in X2). This is for
# Condition 1 (the slower of 2 he tested) and using a Gupta chemistry
# model. We are interested in the VUV region particularly, between
# 100 and 200nm, for comparison with experiment.

# NOTE: because the mass fractions, temperature and pressure are
# specified, we don't need to use CEA, hence the modified EQ_spectra
# script in this folder: EQ_spectra_noCEA.py.

# 1. Setup the radiation model LUA input file
radmodel.py -i air-radiators-photaura.py -L rad-model.lua
# radmodel.py -i air-radiators-parade.py -L rad-model.lua

# 2. Run EQ_spectra.py
EQ_spectra_noCEA.py --input-file=VUV-air.py

# The files coefficient_spectra.txt and intensity_spectra.txt
# should now be in the working directory.
