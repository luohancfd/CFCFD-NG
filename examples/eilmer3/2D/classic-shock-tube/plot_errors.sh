#! /bin/bash
# plot_errors.sh

gnuplot <<EOF
set term postscript eps 20
set output "cst-errors.eps"
set title "Classic shock tube: shock and gas speed with cell resolution."
set xlabel "1/nni"
set ylabel "relative error in U"
set logscale xy
set key right bottom
e(x) = e0 * (nni*x)**n
plot "speeds.data" using (1.0/\$1):(\$2/3603.687-1.0) title "U_s" with points pt 2 ps 1.5, \
     nni = 400, e0 = 0.00960, n = 0.822, e(x) title "n=0.822" lw 2, \
     "speeds.data" using (1.0/\$1):(1.0-\$3/3194.170) title "U_g" with points pt 4 ps 1.5, \
     nni = 400, e0 = 0.00136, n = 0.733, e(x) title "n=0.733" lw 2
EOF