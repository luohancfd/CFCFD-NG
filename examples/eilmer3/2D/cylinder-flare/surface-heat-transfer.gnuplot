# surface-heat-transfer.gnuplot
set term postscript eps 20
set output 'surface-heat-transfer.eps'
set title 'Cylinder with extended flare, heat-flux along surface'
set ylabel 'q, kW/m**2'
set xlabel 'x, mm'
set yrange [0:300]
set key top left
plot './surface.data' using ($1*1000):($9/1000) with lines \
     lw 3.0 title 'Eilmer3 k*dT/dy', \
     './notes/cylinder-extended-flare-heat-transfer.data' \
     using ($2*101.7):($10*11.377) \
     title 'CUBRC Run 14' with points pt 4
