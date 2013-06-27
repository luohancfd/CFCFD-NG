# surface-pressure.gnuplot
set term postscript eps 20
set output 'surface-pressure.eps'
set title 'Double-cone sharp-nose, pressure along surface'
set ylabel 'p, Pa'
set xlabel 'x, mm'
set key top left
plot './surface.data' using ($1*1000):($7) with lines \
     lw 3.0 title 'Eilmer3', \
     './notes/indented-cone-pressure.data' \
     using ($2*92.075):($9*6894.8) \
     title 'CUBRC Run 35' with points pt 4
