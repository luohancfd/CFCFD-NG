#!/bin/sh
# plot-optimised-curves-comparison.sh
#

./generate-optimised-curves-comparison.py

gnuplot<<EOF
set term postscript eps color enhanced size 20cm,20cm
set output "optimised-curves-comparison.eps"
set key right center
set xlabel "Axial distance, m"
set ylabel "Radial distance, m"

#set size ratio -1

plot "optimised-curves-comparison.data" every :::0::0 using (\$1):(\$2) with lines lt 1 lc 2, \
     "" every :::5::5 using (\$1):(\$2) with lines lt 1 lc 6, \
     "" every :::10::10 using (\$1):(\$2) with lines lt 1 lc 3, \
     "" every :::12::12 using (\$1):(\$2) with lines lt 1 lc 4, \
     "" every :::14::14 using (\$1):(\$2) with lines lt 1 lc 5, \
     "" every :::16::16 using (\$1):(\$2) with lines lt 1 lc 1, \
     "../CORRECT_T4_M7_final_20.5mm_fullyTurbulent/optimised-curves-comparison.data" every :::15::15 using (\$1):(\$2) with lines lt 1 lc 7 lw 2, \
     "../CORRECT_T4_M7_final_20.5mm_fullyTurbulent_increase_theta_weighting/optimised-curves-comparison.data" every :::15::15 using (\$1):(\$2) with lines lt 3 lc 7 lw 2
EOF
#     "" every :::15::15 using (\$1):(\$2) with lines lt 1 lc 7 lw 2
#EOF
# plot every :::5::8 plot from 5-th to 8-th data blocks


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

#rm optimised-curves-comparison.data
