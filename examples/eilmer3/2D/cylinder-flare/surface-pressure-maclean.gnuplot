# surface-pressure-maclean.gnuplot
set term postscript eps 20
set output 'surface-pressure-maclean.eps'
set title 'Cylinder with extended flare, pressure along surface'
set ylabel 'p, Pa'
set xlabel 'x, mm'
set key top left
plot './surface.data' using ($1*1000):($7) with lines \
     lw 3.0 title 'Eilmer3', \
     './notes/maclean-2004-0529-figure-6-Cp.data' \
     using ($1*101.7):($2*2338.4) \
     title 'CUBRC Expt' with points pt 4
