# surface-heat-transfer-maclean.gnuplot
set term postscript eps 20
set output 'surface-heat-transfer-maclean.eps'
set title 'Cylinder with extended flare, heat-flux along surface'
set ylabel 'q, kW/m**2'
set xlabel 'x, mm'
set yrange [0:300]
set key top left
plot './surface.data' using ($1*1000):($9/1000) with lines \
     lw 3.0 title 'Eilmer3 k*dT/dy', \
     './notes/maclean-2004-0529-figure-6-Ch.data' \
     using ($1*101.7):($2*5387.6) \
     title 'CUBRC Expt' with points pt 4
