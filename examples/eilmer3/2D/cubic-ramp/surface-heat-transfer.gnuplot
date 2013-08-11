# surface-heat-transfer.gnuplot
set term postscript eps 20
set output 'surface-heat-transfer.eps'
set title 'Cubic ramp, heat-flux along surface'
set ylabel 'q, kW/m**2'
set xlabel 'x, mm'
set yrange [0:200]
set key top left
plot './surface.data' using ($1*1000):($8/1000) with lines \
     lw 3.0 title 'Eilmer3 k*dT/dy', \
     './notes/mohammadian-figure-10-stanton.data' \
     using ($1*25.4):($2*7.682) \
     title 'Mohammadian (1972)' with points pt 4
