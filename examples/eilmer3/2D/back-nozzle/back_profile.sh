# back_profile.sh
# Extract the flow data along the nozzle wall,
# scale it so that it can be lotted with the experimental data
# and plot it using gnuplot.

e3post.py --job=back --output-file=raw_profile.data --tindx=9999 \
    --slice-list=":,:,-1,0"

awk -f normalize.awk raw_profile.data > norm_profile.data

gnuplot <<EOF
set term postscript eps 20
set output 'back_profile_whole.eps'

set title 'Pressure along the nozzle wall'
set xlabel 'distance from throat (inches)'
set ylabel 'p/pt'

set yrange [0:1.2]
plot 'norm_profile.data' using 1:2 title "simulation" with lines, \
     'back-exp.data' using 1:2 title "experiment" with points
EOF

gnuplot <<EOF
set term postscript eps 20
set output 'back_profile_supersonic.eps'

set title 'Pressure along the nozzle wall'
set xlabel 'distance from throat (inches)'
set ylabel 'p/pt'

set xrange [-1.0:3.0]
set yrange [0:0.6]
plot 'norm_profile.data' using 1:2 title "simulation" with lines, \
     'back-exp.data' using 1:2 title "experiment" with points
EOF
