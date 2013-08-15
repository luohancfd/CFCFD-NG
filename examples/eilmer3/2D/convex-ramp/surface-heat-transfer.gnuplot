# surface-heat-transfer.gnuplot
set term postscript eps 20
set output 'surface-heat-transfer.eps'
set title 'Cubic ramp, heat-flux along surface'
set ylabel 'q, kW/m**2'
set xlabel 'x, mm'
set yrange [0:150]
set key top left
plot './surface.data' using ($1*1000):($8/1000) with lines \
     lw 3.0 title 'Eilmer3', \
     './notes/mohammadian-figure-13-heat-flux.data' \
     using ($1*25.4):($2*11.4) \
     title 'Mohammadian (1972)' with points pt 4
