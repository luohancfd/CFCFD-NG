#!/bin/bash
# A test case exercising the new input of atomic radiation data from raw NIST,
# TOPBase and Griem files. 
# A comparison is made between spectra calculated via NIST and TOPBase data.

# 1. NIST data using hydrogenic model for photoionization
radmodel.py -i argon-radiators-NIST.py -L rad-model.lua
EQ_spectra.py --input-file=argon.py
mv coefficient_spectra.txt NIST_coefficient_spectra.txt
mv intensity_spectra.txt NIST_intensity_spectra.txt

# 2. TOPBase data using available tabulated data for photoionization
radmodel.py -i argon-radiators-TOPBase.py -L rad-model.lua
EQ_spectra.py --input-file=argon.py
mv coefficient_spectra.txt TOPBase_coefficient_spectra.txt
mv intensity_spectra.txt TOPBase_intensity_spectra.txt

# 3. Make a plot of the comparison
gnuplot plot.gnu
