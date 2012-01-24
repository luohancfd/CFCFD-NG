#!/bin/sh
# run_simulation.sh

e3prep.py --job=titan_x2_shell 
time e3shared.exe --job=titan_x2_shell --run

echo "Done."
