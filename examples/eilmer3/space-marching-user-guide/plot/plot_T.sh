#! /bin/bash

gnuplot << EOF

set term postscript eps color 20 size 11,8
set output 'T_profile_centreline.eps'

set style line 1 linetype 1 linewidth 7.0 linecolor rgb "black"
set style line 2 linetype 3 linewidth 5.0 linecolor rgb "red"
set style line 3 linetype 3 linewidth 3.0 linecolor rgb "orange"
set style line 4 linetype 3 linewidth 3.0 linecolor rgb "blue"

set title 'Temperature plot down the centreline' font "Helvetica,40"

set key font "Helvetica,30" spacing 1.5

set xlabel 'Distance (m)' font "Helvetica,30"
set ylabel 'Temperature (K)' font "Helvetica,30"

set xrange [0.25:0.85]
set yrange [0.0:2400.0]

plot 'eilmer3_new_code_3_flow_lengths.data' using 1:20 title 'eilmer3; current space marched; 3 flow lengths' with lines linestyle 1,\
'eilmer3_original_3_flow_lengths.data' using 1:20 title 'eilmer3; original space marched; 3 flow lengths' with lines linestyle 4,\
'eilmer3_new_code_6_flow_lengths.data' using 1:20 title 'eilmer3; current space marched; 6 flow lengths' with lines linestyle 3,\
'eilmer3_time_resolved.data' using 1:20 title 'eilmer3; time-dependent' with lines linestyle 2

EOF
