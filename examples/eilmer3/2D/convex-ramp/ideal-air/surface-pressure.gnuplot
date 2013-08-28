# surface-pressure.gnuplot
set term postscript eps 20
set output 'surface-pressure.eps'
set title 'Cubic ramp, pressure along surface'
set ylabel 'p, Pa'
set xlabel 'x, mm'
set key top left
plot './surface.data' using ($1*1000):($6) with lines \
     lw 3.0 title 'Eilmer3', \
     './notes/mohammadian-figure-12-p_p_inf.data' \
     using ($1*25.4):($2*66.43) \
     title 'Mohammadian (1972)' with points pt 4
