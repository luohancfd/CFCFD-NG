#!/bin/sh
# plot-nozzle-profile.sh
#

gnuplot<<EOF
set term postscript eps color enhanced size 20cm,4cm
set output "M7-nozzle-contour.eps"
#set key right center 
unset key
set xlabel "Axial distance, m"
set ylabel "Radial distance, m"
set size ratio -1

plot "../T4_M7_final_20.5mm_fullyTurbulent/finer-grid-run-of-contour-from-iteration-220/contour-t4-m7.data" using (\$1):(\$2) with points pt 1 lc 2, \
     "contour-t4-m7.data" using (\$1):(\$2) with lines lt 1 lc 1
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 20cm,6cm
set output "M7-nozzle-contour-divergence.eps"
#set key right center
unset key
set xlabel "Axial distance, m"
set ylabel "Divergence, degrees"

plot "divergence-t4-m7.data" using (\$1):(\$2) with lines lt 1 lc 1
EOF

epstopdf M7-nozzle-contour.eps
epstopdf M7-nozzle-contour-divergence.eps
rm *eps
