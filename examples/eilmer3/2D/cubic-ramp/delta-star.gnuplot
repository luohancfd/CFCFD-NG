# delta-star.gnuplot
set term postscript eps 20
set output 'delta-star.eps'
set title 'Cubic ramp, edge of boundary layer'
set ylabel 'delta*, mm'
set xlabel 'x, mm'
set key top right
set yrange [0:5]
plot './cubic-ramp-factor-2-grad-rho-field-1-ms-delta.data' \
     using ($1):($2 - 25.4/150.0*($1/25.4)**3) with points pt 1 \
     title 'Eilmer3', \
     './notes/mohammadian-figure-11-delta-star.data' \
     using ($1*25.4):($2*25.4) \
     title 'Mohammadian (1972)' with points pt 4 ps 1.5
