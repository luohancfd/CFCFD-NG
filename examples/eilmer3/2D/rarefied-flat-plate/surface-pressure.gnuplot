# surface-pressure.gnuplot
set term postscript eps 20
set output "surface-pressure.eps"
set title "Rarefied flat plate surface pressure"
set xlabel "x, mm"
set ylabel "static pressure, Pa
set yrange [0:200]
plot "surface.data" using ($1*1000):($9) title "Eilmer3, t=100us" with lines

