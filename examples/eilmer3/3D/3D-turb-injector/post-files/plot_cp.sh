gnuplot<<EOF
set term postscript eps enhanced 20
set output "normalisedCP.eps"
set title "Pressure Coefficient along Centreline"
set xlabel "x/D"
set ylabel "C_P"
set xrange [-10:10]
set yrange [-0.1:0.5]
set key right top
plot "./cp-centreline-complete-31.dat" using (\$2):(\$3) title "Eilmer3 Simulation" with lines lt 1 lw 2, \
     "./centrelineCp-schetz-complete.dat" using (\$1):(\$2) title "Viti Simulation" with lines lt 2 lw 2, \
     "./TC3viti-experimentalCp-complete.dat" using (\$1):(\$2) title "Wallis Data (PSP)" with lines lt 4 lw 2, \
     "./TC3viti-experimentalCp-ERROR.dat" using (\$1):(\$2):(\$3) title "Wallis Error (PSP)" with errorbars
EOF

