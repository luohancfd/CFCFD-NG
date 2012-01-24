#!/bin/bash
# plot.sh
# Plot properties along the stagnation line (and nearby).

gnuplot <<EOF
set term postscript eps enhanced 20
set output "density.eps"
set title "Density along stagnation line, X2 condition, Park-2T."
set xlabel "Normalized position from stagnation point (x/R)"
set ylabel "density, kg/m**3"
set key left bottom
set xrange [0:0.10]
plot "./line0.data" using (-1*\$1/0.02-1.0):(\$5) with lines lt 1 title "line-0", \
     "./line1.data" using (-1*\$1/0.02-1.0):(\$5) with lines lt 2 title "line-1"
EOF