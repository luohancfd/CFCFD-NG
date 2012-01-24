#!/bin/bash

# 1. Single temperature case
poshax.x nitrogen-6.4-singleT.cfg

# 2. Two-temperature case
poshax.x nitrogen-6.4.cfg

# 3. plot results
gnuplot plot_data.gnu
