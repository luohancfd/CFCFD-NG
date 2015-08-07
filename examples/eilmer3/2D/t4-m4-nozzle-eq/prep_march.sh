#! /bin/bash
# prep_march.sh

# Build look-up table for equilibrium air using CEA2.
build-cea-lut.py --gas=air

# Copy up-to-date file containing Bezier control points
# for the T4 Mach 4 nozzle.
cp $HOME/cfcfd3/app/nenzfr/nenzfr_data_files/Bezier-control-pts-t4-m4.data ./

# Run e3prep.
e3prep.py --job=t4-m4-nozzle
