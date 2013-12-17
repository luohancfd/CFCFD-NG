#!/bin/sh
# plot-optimised-curves-comparison.sh
#

./generate-initial-contour.py
./generate-optimised-curves-comparison.py

gnuplot<<EOF
set term postscript eps color enhanced size 20cm,4cm
set output "optimised-curves-comparison.eps"
set key right center
set xlabel "Axial distance, m"
set ylabel "Radial distance, m"
set size ratio -1

plot "contour-t4-m4.initial.data" using (\$1):(\$2) title "Initial contour" with lines lt 1 lc 1, \
     "optimised-curves-comparison.data" every :::0::0 using (\$1):(\$2) title "after 0-th optimisation run" with lines lt 1 lc 2, \
     "" every :::5::5 using (\$1):(\$2) title "after 5-th optimisation run" with lines lt 1 lc 3, \
     "" every :::10::10 using (\$1):(\$2) title "after 10-th optimisation run" with lines lt 1 lc 4, \
     "" every :::12::12 using (\$1):(\$2) title "after 12-th optimisation run" with lines lt 1 lc 5, \
     "" every :::14::14 using (\$1):(\$2) title "after 14-th optimisation run" with lines lt 1 lc 6
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "residuals.eps"
set logscale y
set xlabel "Number of iterations"
set ylabel "Objective function value"
unset key
plot "residuals.data" using (\$1):(\$3) with linespoint lt 1 lc 1
EOF

epstopdf optimised-curves-comparison.eps
epstopdf residuals.eps
rm *eps
