#!/bin/bash
# post_simulation.sh

# Extract the stagnation line data from the steady flow field.
e3post.py --job=bar --output-file=stag_line.data --tindx=9999 \
     --slice-list="0,:,1,0"
awk -f locate_shock.awk stag_line.data > result.txt

# Create a VTK plot file of the steady flow field.
e3post.py --job=bar --tindx=all --vtk-xml --add-mach --add-pitot-p

# Extract the flow data across the face of the bar gauge.
e3post.py --job=bar --output-file=raw_profile.data --tindx=9999 \
     --slice-list="0,-1,:,0"
awk -f normalize.awk raw_profile.data > norm_profile.data

gnuplot <<EOF
set output "bar_norm_p.eps"
set term postscript eps 20
set xrange [0:1.1]
set yrange [0:1.2]
set title "Normalized surface pressure over cylinder face, M=4.76."
set xlabel "r/r-max"
set ylabel "p/p-centre"
set key bottom left
plot "norm_profile.data" using 1:2 title "simulation" with lines, \
     "kendall_profile.data" using 1:2 title "experiment" with points pt 4
EOF
