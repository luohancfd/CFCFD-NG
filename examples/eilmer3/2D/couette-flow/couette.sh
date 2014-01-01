#!/bin/sh
# couette.sh

e3prep.py --job=couette --do-svg
e3post.py --job=couette --vtk-xml --tindx=0

time e3shared.exe --job=couette --run

e3post.py --job=couette --vtk-xml --tindx=last

e3post.py --job=couette --output-file=dudy0.dat --tindx=0 \
           --slice-list="0,1,:,0"
e3post.py --job=couette --output-file=dudy1.dat --tindx=last \
           --slice-list="0,1,:,0"

gnuplot <<EOF
set term postscript eps 20
set output "velocity.ps"
set title "Velocity profile along the height"
set ylabel "Height, m"
set xlabel "Velocity, m/s"
set yrange [0.0:0.0115]
set xrange [-10.0:110.0]
plot "dudy0.dat" using 6:2 with lines title "Initial value", \
     "dudy1.dat" using 6:2 with lines title "Steady state condition"
EOF

