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
plot "speeds.data" using (1.0/\$1):(\$2/3603.687-1.0) title "U_s" with points pt 2, \
     "speeds.data" using (1.0/\$1):(1.0-\$3/3194.170) title "U_g" with points pt 4
EOF