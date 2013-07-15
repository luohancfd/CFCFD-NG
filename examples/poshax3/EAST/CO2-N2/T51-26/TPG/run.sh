#!/bin/bash

# 1. Prepare gas and radiation files
gasfile mars-TPG.inp gas-model.lua
radmodel.py -i mars-radiators.py -L rad-model.lua

# 2. Run simulation
poshax3.x EAST.cfg

# 3. Plot results
gnuplot profiles.gplot
gnuplot spectra.gplot
