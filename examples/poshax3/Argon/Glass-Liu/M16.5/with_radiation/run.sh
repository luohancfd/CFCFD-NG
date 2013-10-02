#!/bin/bash

# 1. Prepare gas and radiation files
gasfile argon-3sp-2T.inp argon-3sp-2T.lua
radmodel.py -i argon-radiators-NIST-TB.py -L rad-model.lua

# 2. Run simulation
poshax3.x argon.cfg

# 3. Plot results
gnuplot profiles.gplot

# 4. Convert to png files
convert -density 600x600 -quality 90 temperature_profiles.eps -scale 750 temperature_profiles.png 
convert -density 600x600 -quality 90 ionization_fraction_profile.eps -scale 750 ionization_fraction_profile.png

# 5. Compute RMS error with UTIAS shock tube measurements
python compute_errors.py 
