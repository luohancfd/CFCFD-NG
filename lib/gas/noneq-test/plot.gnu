set term postscript eps enhanced color
set grid
set size 1.5,1.0
set key top left

set output "temperature_compare.eps"
set xlabel "Evaluation number"
set ylabel "Temperature, T (K)"
plot '5sp-EQ_libgas_results.txt' u 1:2 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 1:2 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 1:2 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 1:2 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 1:2 w l lw 3 t 'cea2 - 11sp'

set xlabel "Temperature, T (K)"

set output "density_compare.eps"
set ylabel "Density, rho (kg/m**3)"
plot '5sp-EQ_libgas_results.txt' u 2:3 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 2:3 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 2:3 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 2:3 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 2:3 w l lw 3 t 'cea2 - 11sp'
     
set output "enthalpy_compare.eps"
set ylabel "Enthalpy, h (kJ/kg)"
plot '5sp-EQ_libgas_results.txt' u 2:4 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 2:4 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 2:4 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 2:4 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 2:4 w l lw 3 t 'cea2 - 11sp'
     
set output "energy_compare.eps"
set ylabel "Energy, e (kJ/kg)"
plot '5sp-EQ_libgas_results.txt' u 2:5 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 2:5 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 2:5 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 2:5 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 2:5 w l lw 3 t 'cea2 - 11sp'
     
set output "entropy_compare.eps"
set ylabel "Entropy, s (kJ/kg-K)"
plot '5sp-EQ_libgas_results.txt' u 2:6 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 2:6 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 2:6 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 2:6 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 2:6 w l lw 3 t 'cea2 - 11sp'
     
set output "Cp_compare.eps"
set ylabel "Cp (kJ/kg-K)"
plot '5sp-EQ_libgas_results.txt' u 2:7 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 2:7 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 2:7 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 2:7 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 2:7 w l lw 3 t 'cea2 - 11sp'
     
set output "gamma_compare.eps"
set ylabel "Specific heat ratio, gamma (ND)"
plot '5sp-EQ_libgas_results.txt' u 2:8 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 2:8 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 2:8 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 2:8 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 2:8 w l lw 3 t 'cea2 - 11sp'
     
set output "son_compare.eps"
set ylabel "Sound speed, a (m/s)"
plot '5sp-EQ_libgas_results.txt' u 2:9 w l lw 3 t 'lib/gas - 5sp-EQ', \
     '5sp-NEQ_libgas_results.txt' u 2:9 w l lw 3 t 'lib/gas - 5sp-NEQ', \
     '5sp-EQ_cea2_results.txt'   u 2:9 w l lw 3 t 'cea2 - 5sp', \
     '11sp-NEQ_libgas_results.txt' u 2:9 w l lw 3 t 'lib/gas - 11sp-NEQ', \
     '11sp-NEQ_cea2_results.txt'   u 2:9 w l lw 3 t 'cea2 - 11sp'
     