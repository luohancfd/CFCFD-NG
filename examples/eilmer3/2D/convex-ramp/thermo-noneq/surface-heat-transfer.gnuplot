# surface-heat-transfer.gnuplot
set term postscript eps 20
set output 'surface-heat-transfer.eps'
set title 'Cubic ramp, heat-flux along surface'
set ylabel 'q, kW/m**2'
set xlabel 'x, mm'
set yrange [0:150]
set key top left
plot './my-surface.data' using ($10*1000):($2/1000) with lines \
     lw 3.0 title 'Eilmer3', \
     './my-surface.data' using ($10*1000):($2/1000*1.2) with lines \
     lw 3.0 lt 2 title 'Eilmer3*1.2', \
     './notes/mohammadian-figure-13-heat-flux.data' \
     using ($1*25.4):($2*11.4) \
     title 'Mohammadian (1972) expt' with points pt 4, \
     './notes/mohammadian-figure-13-original-cheng-theory.data' \
     using ($1*25.4):($2*11.4) \
     title 'original Cheng theory' with lines lw 1.5 lt 3, \
     './notes/mohammadian-figure-13-modified-cheng-theory.data' \
     using ($1*25.4):($2*11.4) \
     title 'modified Cheng theory' with lines lw 1.5 lt 4
