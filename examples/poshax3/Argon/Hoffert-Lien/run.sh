#!/bin/bash

# 1. Prepare gas file
gasfile argon-3sp-2T.inp argon-3sp-2T.lua

# 2. Run simulation
poshax3.x argon.cfg

# 3. Plot results
gnuplot profiles.gplot

# 4. Convert to png files
# convert -density 600x600 -quality 90 number_density_profiles.eps number_density_profiles.png
# convert -density 600x600 -quality 90 temperature_profiles.eps temperature_profiles.png 

# 5. Compute average errors with Panesi solution
# python compute_errors.py 
