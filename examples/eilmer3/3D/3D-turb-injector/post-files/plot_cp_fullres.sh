gnuplot<<EOF
set term postscript eps enhanced 20
set output "normalisedCP_fullres.eps"
set title "Pressure Coefficient along Centreline"
set xlabel "x/D"
set ylabel "C_P"
set xrange [-16:44]
set yrange [-0.1:0.5]
set key right top
plot "./cp-centreline-complete-31.dat" using (\$2):(\$3) title "Eilmer3 Simulation" with lines lt 1 lw 2, \
     "./centrelineCp-schetz-complete.dat" using (\$1):(\$2) title "Viti Simulation" with lines lt 2 lw 2
EOF

