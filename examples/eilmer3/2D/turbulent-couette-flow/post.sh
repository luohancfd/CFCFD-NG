#!/bin/sh
e3post.py --job=couette --vtk-xml --tindx=last

e3post.py --job=couette --bc-surface-list="0,3,0,:,0" --tindx=last

e3post.py --job=couette --output-file="slice.txt" --slice-list="42,0,:,0" --tindx=last

gnuplot <<EOF
set grid
set term postscript eps 20
set output "boundary-layer.eps"
set xlabel "Non-dimensional height"
set ylabel "Non-dimensional velocity"
set xrange [0.0:2.0]
set yrange [0.0:1.0]
plot "slice.txt" using (2.0-\$2/33.0e-3):(\$6/12.84) with lines title "Eilmer3", \
     "telbany.dat" using 1:2 with points title "Telbany's experiment"
EOF

gnuplot <<EOF
set key left
set grid x y mx my
set logscale x
set term postscript eps 20
set output "law-of-the-wall.eps"
set xlabel "y+"
set ylabel "u+"
set xrange [0.1:10000.0]
set yrange [0.0:28.0]
f(x) = x
y(x) = 1.0/0.41*log(x)+5.5
plot f(x) title "Viscous sublayer", \
     y(x) title "log-law region", \
     "eilmer.txt" using 1:2 with linespoints title "Eilmer3"
EOF
