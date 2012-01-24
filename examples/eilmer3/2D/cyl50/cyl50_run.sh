#! /bin/sh
# cyl50_run.sh
# exercise the Navier-Stokes solver for the cyl50 test case.
# It is assumed that the path is set correctly.

# Generate the Input parameter, grid, initial flow and SVG files.
# The SVG file provides us with a graphical check on the geometry
# specification and it is worth creating if a viewer is available.
e3prep.py --job=cyl50 --do-svg

# Integrate the solution in time, 
time e3shared.exe --job=cyl50 --run

echo "At this point, we should have a new solution"
echo "Run cyl50_plot.sh next"

