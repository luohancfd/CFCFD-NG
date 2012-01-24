# back_fig2.gnu
set output 'back_ppt.ps'
set term postscript eps
set xlabel 'x, inches'
set ylabel 'p/pt'
set xrange [0:3.0]
set yrange [0:0.6]
plot 'back_exp.dat' using 1:2, 'back_moc.dat' using 1:2
