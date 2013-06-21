#!/bin/bash
# This script plots the results from the Photaura simulation.

# 1. Plot the results
gnuplot plot-spectra.gnu

# 2. Convert plots to .pdf files and view plots
epstopdf Intensity-spectra.eps
evince Intensity-spectra.pdf

# The plots should be in the working directory.
