#!/bin/bash

# 1. Setup the radiation model LUA input file
radmodel.py -i H2-air-radiators.py -L rad-model.lua

# 2. Run EQ_spectra.py
EQ_spectra.py --input-file=H2-air.py

# The files coefficient_spectra.txt and intensity_spectra.txt
# should now be in the working directory.
