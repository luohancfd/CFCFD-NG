#!/bin/sh
# plot-optimised-curves-comparison.sh
#

./generate-initial-contour.py

gnuplot<<EOF
set term postscript eps color enhanced size 20cm,20cm
set output "compare-initial-contours.eps"
set key right center
set xlabel "Axial distance, m"
set ylabel "Radial distance, m"

#set size ratio -1

plot "contour-t4-m7.initial.data" using (\$1):(\$2) title "Contour scaled from M8 nozzle contour" with lines lt 1 lc 1 lw 2, \
     "../T4_M8_test/contour-t4-m8.data" using (\$1):(\$2) title "M8 nozzle contour" with lines lt 2 lc 3 lw 2, \
     "../T4_M7_final_20.5mm_fullyTurbulent/optimised-curves-comparison.data" every :::0::0 using (\$1):(\$2) title "IMOC scaled contour" with lines lt 5 lc 1 lw 2
EOF

epstopdf compare-initial-contours.eps
rm *eps
