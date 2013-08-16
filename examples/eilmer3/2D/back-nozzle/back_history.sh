# back_history.sh
# Extract the flow history data at the nozzle exit plane.
# This is then plotted using gnuplot and an assessment
# can be made as to whether the flow has reached steady state.

awk -f extract-history.awk < hist/back.hist.b0001 > nozzle-exit.data

gnuplot <<EOF
set term postscript eps 20
set output 'back_history_M.eps'

set title 'Mach number history at the nozzle exit'
set xrange [0.0:4.0]
set xlabel 'time, ms'
set ylabel 'M'

plot 'nozzle-exit.data' using 1:2 with lines
EOF

gnuplot <<EOF
set term postscript eps 20
set output 'back_history_p.eps'

set title 'Static pressure history at the nozzle exit'
set key bottom right
set xrange [0.0:4.0]
set xlabel 'time, ms'
set ylabel 'p, kPa'

plot 'nozzle-exit.data' using 1:3 with lines
EOF

