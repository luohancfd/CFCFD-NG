# surface-pressure.gnuplot
set term postscript eps 20
set output 'surface-pressure.eps'
set title 'Cylinder with extended flare, pressure along surface'
set ylabel 'p, Pa'
set xlabel 'x, mm'
set key top left
plot './surface.data' using ($1*1000):($7) with lines \
     lw 3.0 title 'Eilmer3', \
     './notes/cylinder-extended-flare-pressure.data' \
     using ($2*101.7):($10*6894.8) \
     title 'CUBRC Run 14' with points pt 4
