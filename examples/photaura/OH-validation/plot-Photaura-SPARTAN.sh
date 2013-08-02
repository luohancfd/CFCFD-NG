#!/bin/bash
# This script plots the results from the Photaura simulation.

# 1. Plot the results
gnuplot plot-Photaura-SPARTAN.gnu

# 2. Convert plots to .pdf files and view plots
epstopdf Intensity-Photaura-SPARTAN.eps
evince Intensity-Photaura-SPARTAN.pdf

# The plots should be in the working directory.
