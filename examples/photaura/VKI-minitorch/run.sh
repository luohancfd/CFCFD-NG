#!/bin/bash

# 1. Prepare the gas model
gasfile N2-CH4-gas.inp gas-model.lua

# 2. Prepare the radiation model
script_rad2.py -i TC4-radiators.py -L TC4-radiators.lua

# 3. Run the calculation script
./TC4-validation.py

# 4. Plot the results
gnuplot plot_spectra.gnu

# 5. Convert to png images for sphinx documentation
for X in Emissivity-r1.24mm Emissivity-r3.82mm Emissivity-r8.43mm Intensity3500-4300A Intensity4300-10000A
do
convert $X.eps $X.png
done
