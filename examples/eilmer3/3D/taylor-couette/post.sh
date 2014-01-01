#!/bin/sh
# post-processing script
# post.sh

gnuplot <<EOF
set term postscript eps 20
set output "velocity.eps"
set title "axially averaged tangential velocity profile"
set xlabel "radial position"
set ylabel "velocity"
set xrange [0.0:1.0]
set yrange [0.0:1.0]
plot "average.txt" using 1:2 with lines title "Eilmer3", \
     "tangential.dat" using 1:2 with lines title "CTDNS"
EOF

gnuplot <<EOF
set term postscript eps 20
set output "temperature.eps"
set title "axially averaged temperature profile"
set xlabel "radial position"
set ylabel "temperature"
set xrange [0.0:1.0]
set yrange [350.0:400.0]
plot "average.txt" using 1:3 with lines title "Eilmer3", \
     "temperature.dat" using 1:2 with lines title "CTDNS"
EOF
