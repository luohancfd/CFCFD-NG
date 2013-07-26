#!/bin/bash

# 1. Prepare gas and radiation files
gasmodel.py --model='two temperature gas' --species='H2 H H_plus He e_minus' --output='gas-model.lua'
radmodel.py -i H2-He-radiators.py -L rad-model.lua

# 2. Run simulation
poshax3.x x2_H2-He.cfg

# 3. Plot results
gnuplot profiles.gplot
gnuplot spectra.gplot
