gnuplot<<EOF
set term postscript eps enhanced 20
set output "normalisedCP_all.eps"
set title "Pressure Coefficient along Centreline"
set xlabel "x/D"
set ylabel "C_P"
#set xrange [-10:10]
set yrange [-0.1:0.5]
set key right top
plot "./cp-centreline-complete-28.dat" using (\$2):(\$3) title "Eilmer3 t=6.828827e-04" with lines lt 4 lw 2, \
     "./cp-centreline-complete-29.dat" using (\$2):(\$3) title "Eilmer3 t=7.078827e-04" with lines lt 3 lw 2, \
     "./cp-centreline-complete-30.dat" using (\$2):(\$3) title "Eilmer3 t=7.328828e-04" with lines lt 2 lw 2, \
     "./cp-centreline-complete-31.dat" using (\$2):(\$3) title "Eilmer3 t=7.578828e-04" with lines lt 1 lw 2
     
EOF

