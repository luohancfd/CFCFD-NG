#!/bin/bash

# 1. Prepare gas file
gasfile air-11sp-2T.inp air-11sp-2T.lua
# We want harmonic oscillators instead of truncated harmonic
# oscillators for this literature comparison.
sed -ie 's/truncated harmonic/harmonic/g' air-11sp-2T.lua
mv air-11sp-2T.lua air-11sp-2T-HO.lua

# 2. Run simulation
poshax3.x FireII.cfg

# 3. Plot results
gnuplot profiles.gplot

# 4. Convert to png files
convert number_density_profiles.eps number_density_profiles.png
convert temperature_profiles.eps temperature_profiles.png 

# 5. Compute average errors with Panesi solution
python compute_errors.py 
