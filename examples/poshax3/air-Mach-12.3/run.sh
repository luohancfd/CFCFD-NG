#!/bin/bash

# 1. Prepare gas file
gasfile air-model.inp air-model.lua

# 2. Run simulation
poshax3.x air-Mach-12.3.cfg

# 3. Plot results
gnuplot plot-data.gplot
convert profile_moles.eps profile_moles.png
convert profile_T_rho.eps profile_T_rho.png

