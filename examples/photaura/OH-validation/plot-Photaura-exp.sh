#!/bin/bash
# This script plots the results from the Photaura simulation.

# 1. Plot the results
gnuplot plot-Photaura-exp.gnu

# 2. Convert plots to .pdf files and view plots
epstopdf Intensity-spectra-Photaura-exp.eps
evince Intensity-spectra-Photaura-exp.pdf

# The plots should be in the working directory.
