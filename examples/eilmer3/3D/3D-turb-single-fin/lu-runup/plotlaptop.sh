gnuplot<<EOF
set term postscript eps enhanced 20
set output "lu-cf-runup.eps"
set title "Lu Flat Plate - c_f profile along wall"
set xlabel "x, m"
set ylabel "c_f"
set yrange [0:0.0025]
set key right top
plot "./viscous_data.dat" using (\$1):(\$5) title "Eilmer3 Lu Input" with lines lt 1 lw 2, \
     "./experimental_data_cf.txt" using (\$1):(\$2) \
     title "Lu Data" with points pt 6 lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20    
set output "lu-mach-runup.eps"
set title "Lu Flat Plate - Mach number profile near wall"
set xlabel "Mach number"
set ylabel "y, mm"
set key left top Left
set xrange [0:5]
set yrange [0:5] #25
plot "./flat_runup-178mm-slice.dat" using (\$21):((\$2)*1000) title "Eilmer3 Lu Input at 0.178m" with lines lt 1 lw 2, \
"./flat_runup-178mm-slice3.dat" using (\$21):((\$2)*1000) title "Eilmer3 Lu Input at 0.197m" with lines lt 2 lw 1, \
     "./experimental_data.txt" using (\$3):((\$2)*1000) \
     title "Lu Data" with point pt 6 lt 1
EOF
