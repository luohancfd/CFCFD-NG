#!/bin/sh
# plot-exit-profiles.sh
#

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "profile-Mach.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Mach number"
set xrange [0:0.135]
set yrange [0:8]
plot "nozzle-exit.data" using (\$2):(\$21) title "from-scaled-M8-nozzle" with lines lt 1 lc 1, \
     "../T4_M7_final_20.5mm_fullyTurbulent/fine-grid-run-of-contour-from-iteration-220/nozzle-exit.data" using (\$2):(\$21) title "from-IMOC-contour" with lines lt 3 lc 1, \
     "../CORRECT_T4_M7_final_20.5mm_fullyTurbulent_start_from_scaled_M8_nozzle_increase_throatDia/nozzle-exit.data" using (\$2):(\$21) title "higher-theta-weighting" with lines lt 5 lc 1
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "profile-Pitot.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Pitot pressure, kPa"
set xrange [0:0.135]
set yrange [0:100]
plot "nozzle-exit.data" using (\$2):((\$22)/1000) title "from-scaled-M8-nozzle" with lines lt 1 lc 1, \
     "../T4_M7_final_20.5mm_fullyTurbulent/fine-grid-run-of-contour-from-iteration-220/nozzle-exit.data" using (\$2):((\$22)/1000) title "from-IMOC-contour" with lines lt 3 lc 1, \
     "../CORRECT_T4_M7_final_20.5mm_fullyTurbulent_start_from_scaled_M8_nozzle_increase_throatDia/nozzle-exit.data" using (\$2):((\$22)/1000) title "higher-theta-weighting" with lines lt 5 lc 1
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "profile-flow-divergence.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Flow divergence, degrees"
set xrange [0:0.135]
set yrange [-5:5]
plot "nozzle-exit.data" using (\$2):(atan((\$7)/(\$6)*180/pi)) title "from-scaled-M8-nozzle" with lines lt 1 lc 1, \
     "../T4_M7_final_20.5mm_fullyTurbulent/fine-grid-run-of-contour-from-iteration-220/nozzle-exit.data" using (\$2):(atan((\$7)/(\$6)*180/pi)) title "from-IMOC-contour" with lines lt 3 lc 1, \
     "../CORRECT_T4_M7_final_20.5mm_fullyTurbulent_start_from_scaled_M8_nozzle_increase_throatDia/nozzle-exit.data" using (\$2):(atan((\$7)/(\$6)*180/pi)) title "higher-theta-weighting" with lines lt 5 lc 1
EOF

epstopdf profile-Mach.eps
epstopdf profile-Pitot.eps
epstopdf profile-flow-divergence.eps
rm *eps
