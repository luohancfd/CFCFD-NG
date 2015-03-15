#!/bin/bash

# 1. Preparation
e3prep.py --job=reactor --clean-start
rm -rf hist/

# 2. Simulation
e3shared.exe --job=reactor --run

# 3. Post-processing
e3history.py --cell=1 --output=cpu-chem-hist.dat hist/reactor.hist.b0000
gnuplot plot-history.gplot


