#! /bin/bash
# make-plots.sh

# Plot the OH and H2O plots as they develop down the tube 
# and compare with the reported Bittker-Scullin data.

gnuplot << EOF
set term postscript eps enhanced 20
set output 'OH_production.eps'
set style line 1 linewidth 3.0
set style line 16 pointsize 2.0
set title 'OH mass fraction along duct'
set key top left
set xlabel 'Distance (m)'
set ylabel 'Mass fraction'
set xrange [0.0:0.1]
set yrange [0.0:0.04]
plot 'eilmer3-computed.data' using 1:25 title 'Eilmer3' with lines linestyle 1,\
     'reacting-pipe-flow.transcript' using 1:7 title 'pipe-flow' with lines linestyle 2, \
     'bittker-scullin-massf.data' using 1:3 title 'Bittker-Scullin' with points linestyle 16;
EOF

gnuplot << EOF
set term postscript eps enhanced 20
set output 'H2O_production.eps'
set style line 1 linewidth 3.0
set style line 16 pointsize 2.0
set title 'H_2O mass-fraction along duct'
set key top left
set xlabel 'Distance (m)'
set ylabel 'Mass fraction'
set xrange [0.0:0.1]
set yrange [0.0:0.20]
plot 'eilmer3-computed.data' using 1:23 title 'Eilmer3' with lines linestyle 1,\
     'reacting-pipe-flow.transcript' using 1:8 title 'pipe-flow' with lines linestyle 2, \
     'bittker-scullin-massf.data' using 1:4 title 'Bittker-Scullin' with points linestyle 16;
EOF

gnuplot << EOF
set term postscript eps enhanced 20
set output 'T-development.eps'
set style line 1 linewidth 3.0
set style line 16 pointsize 2.0
set title 'Static temperature along duct'
set key top left
set xlabel 'Distance (m)'
set ylabel 'T (K)'
set xrange [0.0:0.1]
set yrange [1500.0:3000.0]
plot 'eilmer3-computed.data' using 1:29 title 'Eilmer3' with lines linestyle 1,\
     'reacting-pipe-flow.transcript' using 1:4 title 'pipe-flow' with lines linestyle 2, \
     'bittker-scullin-massf.data' using 1:11 title 'Bittker-Scullin' with points linestyle 16;
EOF

gnuplot << EOF
set term postscript eps enhanced 20
set output 'p-development.eps'
set style line 1 linewidth 3.0
set style line 16 pointsize 2.0
set title 'Static pressure along duct'
set key top left
set xlabel 'Distance (m)'
set ylabel 'p (kPa)'
set xrange [0.0:0.1]
set yrange [90.0:160.0]
plot 'eilmer3-computed.data' using (\$1):(\$9/1000) title 'Eilmer3' with lines linestyle 1,\
     'reacting-pipe-flow.transcript' using (\$1):(\$3/1000) title 'pipe-flow' with lines linestyle 2, \
     'bittker-scullin-massf.data' using (\$1):(\$12/1000) title 'Bittker-Scullin' with points linestyle 16;
EOF
