#!/bin/sh
# couette.sh

e3prep.py --job=couette --do-svg
e3post.py --job=couette --vtk-xml --tindx=0

time e3shared.exe --job=couette --run

e3post.py --job=couette --vtk-xml --tindx=last

# post-processing
e3post.py --job=couette --output-file=dudz0.dat --tindx=0 \
           --slice-list="0,-1,10,:"
e3post.py --job=couette --output-file=dudz5.dat --tindx=5 \
           --slice-list="0,-1,10,:"

gnuplot <<EOF
set term postscript eps 20
set output "Velocity.ps"
set title "Velocity profile along the height"
set ylabel "Height, m"
set xlabel "Velocity, m/s"
set yrange [0.0:3.5e-3]
set xrange [-10.0:140.0]
plot "dudz0.dat" using 7:3 with lines title "Initial value", \
     "dudz5.dat" using 7:3 with lines title "Steady state condition"
EOF
