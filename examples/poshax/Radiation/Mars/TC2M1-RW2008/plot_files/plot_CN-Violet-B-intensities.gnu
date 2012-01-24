set term postscript eps enhanced color "Helvetica" 28
set size 1.75, 1.5
set xlabel "Distance behind shock, x (mm)"
set ylabel "Intensity, I_{/Symbol l} (W/cm^2-sr)"
set output "CN-Violet-B_IvX.eps"
set xrange [-2:21.0]
set grid
plot "2T-basic-Lee-CN-Violet-B_IvX.dat" u ($1*1000.0):($2*1.0e-4)  t '2T Basic Lee' with linespoints lt 2 lw 8, \
     "2T-full-Lee-CN-Violet-B_IvX.dat"  u ($1*1000.0):($2*1.0e-4)  t '2T Full Lee' with linespoints lt 5 lw 8, \
     "2T-full-Park-CN-Violet-B_IvX.dat"  u ($1*1000.0):($2*1.0e-4) t '2T Full Park' with linespoints lt 8 lw 8, \
     "2T-full-PL-CN-Violet-B_IvX.dat"   u ($1*1000.0):($2*1.0e-4)  t '2T Full PL' with linespoints lt 4 lw 8, \
     "2T-full-PLF-CN-Violet-B_IvX.dat"  u ($1*1000.0):($2*1.0e-4)  t '2T Full PLF' with linespoints lt -1 lw 4, \
     "3T-full-PLF-CN-Violet-B_IvX.dat"  u ($1*1000.0):($2*1.0e-4)  t '3T Full PLF' with linespoints lt 0 lw 8, \
     "2T-full-PLFM-CN-Violet-B_IvX.dat" u ($1*1000.0):($2*1.0e-4)  t '2T Full PLFM' with linespoints lt 1 lw 8, \
     "2T-full-PLF-RC-NEQ-CN-Violet-B_IvX.dat" u ($1*1000.0):($2*1.0e-4)   t '2T Full PLF - RC NEQ' with linespoints lt 3 lw 8, \
     "2T-full-PLF-RC-EQ-CN-Violet-B_IvX.dat" u ($1*1000.0):($2*1.0e-4)   t '2T Full PLF - RC EQ' with linespoints lt 7 lw 8, \
     "EAST-CN-Violet-B_IvX.txt" u ($1*-10+76):2 t "EAST" w p lt -1 ps 4 pt 1
