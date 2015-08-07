#!/bin/sh
# plot-exit-profiles.sh
#

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "profile-Mach.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Mach number"
set xrange [0:0.08]
#set yrange [0:5]
set key bottom left
plot "nozzle-exit.data" using (\$2):(\$21) title "e3march-version" with lines lt 1 lc 1 lw 2, \
     "nozzle-exit.nenzfr.data" using (\$2):(\$21) title "NENZFr-version" with lines lt 1 lc 3 lw 2
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "profile-Pitot.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Pitot pressure, kPa"
set xrange [0:0.08]
#set yrange [0:700]
set key bottom left
plot "nozzle-exit.data" using (\$2):((\$22)/1000) title "e3march-version" with lines lt 1 lc 1 lw 2, \
     "nozzle-exit.nenzfr.data" using (\$2):((\$22)/1000) title "NENZFr-version" with lines lt 1 lc 3 lw 2
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "profile-staticP.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Static pressure, kPa"
set xrange [0:0.08]
#set yrange [0:20]
set key bottom left
plot "nozzle-exit.data" using (\$2):((\$9)/1000) title "e3march-version" with lines lt 1 lc 1 lw 2, \
     "nozzle-exit.nenzfr.data" using (\$2):((\$9)/1000) title "NENZFr-version" with lines lt 1 lc 3 lw 2
EOF

gnuplot<<EOF
set term postscript eps color enhanced size 12cm,9cm
set output "profile-flow-divergence.eps"
set xlabel "Radial distance from axis, m"
set ylabel "Flow divergence, degrees"
set xrange [0:0.08]
set yrange [-2:2]
set key bottom left
plot "nozzle-exit.data" using (\$2):(atan((\$7)/(\$6)*180/pi)) title "e3march-version" with lines lt 1 lc 1 lw 2, \
     "nozzle-exit.nenzfr.data" using (\$2):(atan((\$7)/(\$6)*180/pi)) title "NENZFr-version" with lines lt 1 lc 3 lw 2
EOF

epstopdf profile-Mach.eps
epstopdf profile-Pitot.eps
epstopdf profile-staticP.eps
epstopdf profile-flow-divergence.eps
