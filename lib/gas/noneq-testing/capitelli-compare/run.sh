#!/bin/bash

for X in N2 NO O2 CO2 CO C N O
do
./noneq-thermo-test.py $X

gnuplot <<EOF
set term postscript eps enhanced "Helvetica,16"
set grid
set size 1.5,1.0
set key top left

set xlabel "Temperature, T (K)"
     
set output "Cp_compare.eps"
set ylabel "Specific heat at constant pressure, C_p (J/mol-K)"
plot 'capitelli_data/$X.txt' u 1:2 w l lw 3 t 'Capitelli (2006)', \
     'libgas_results/$X.txt'  u 1:2 w l lw 3 t 'lib/gas'
     
set output "Cv_compare.eps"
set ylabel "Specific heat at constant volume, C_v (J/mol-K)"
plot 'capitelli_data/$X.txt' u 1:4 w l lw 3 t 'Capitelli (2006)', \
     'libgas_results/$X.txt'  u 1:3 w l lw 3 t 'lib/gas'
     
set output "gamma_compare.eps"
set ylabel "Specific heat ratio, {/Symbol g} (ND)"
plot 'capitelli_data/$X.txt' u 1:5 w l lw 3 t 'Capitelli (2006)', \
     'libgas_results/$X.txt'  u 1:4 w l lw 3 t 'lib/gas'
EOF

mkdir plots/$X

mv *.eps plots/$X/
done
