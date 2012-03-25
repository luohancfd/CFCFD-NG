#!/bin/bash

# Extract the profile along the shock tube.
# slice-list=block-range,i-range,j-range,k-range
#            all blocks - :
#            all i's - :
#            constant j - 0
#            constant k - 0 (not relevant in 2D anyway)
e3post.py --job=cst --slice-list=":,:,0,0" --output-file=profile.data

# Plot the data along the x-axis.
gnuplot <<EOF
set term postscript eps enhanced 20
set output "cst-p.eps"
set title "High-performance shock tube at t = 100us"
set xlabel "x, m"
set ylabel "Pressure, MPa"
set xrange [0:1]
set yrange [0:32.0]
plot "profile.data" using 1:(\$9/1.0e6) t 'Eilmer3' with points pt 6, \
     "exact.data" using 1:(\$3/1.0e6) t 'analytic soln' with lines lt 1 lw 3
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "cst-rho.eps"
set title "High-performance shock tube at t = 100us"
set xlabel "x, m"
set ylabel "Density, kg/m^3"
set xrange [0:1]
set yrange [0:5]
set key left bottom
plot "profile.data" using 1:5 t 'Eilmer3' with points pt 6, \
     "exact.data" using 1:2 t 'analytic soln' with lines lt 1 lw 3
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "cst-u.eps"
set title "High-performance shock tube at t = 100us"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set xrange [0:1]
set yrange [0:3500]
set key left top
plot "profile.data" using 1:6 t 'Eilmer3' with points pt 6, \
     "exact.data" using 1:5 t 'analytic soln' with lines lt 1 lw 3
EOF

gnuplot <<EOF
set term postscript eps enhanced 20
set output "cst-T.eps"
set title "High-performance shock tube at t = 100us"
set xlabel "x, m"
set ylabel "Temperature, degrees K"
set xrange [0:1]
set yrange [0:5000]
set key left bottom
plot "profile.data" using 1:22 t 'Eilmer3' with points pt 6, \
     "exact.data" using 1:4 t 'analytic soln' with lines lt 1 lw 3
EOF



