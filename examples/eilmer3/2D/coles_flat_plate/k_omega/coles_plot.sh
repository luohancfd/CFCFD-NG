#!/bin/sh
# coles_plot.sh
#
# Plot the profile data at approx. x=546mm so that we can compare
# the simulation results with data from Coles 1953 experiment.
# The dimensional data from the experiment is reconstructed from 
# the normalised data recorded in the AGARD report 223. 
#
# Peter J., 30-Oct-2007 
#           31-Aug-2008 updated for Elmer3
#

gnuplot<<EOF
set term postscript eps 20
set output "coles-x05-u.eps"
set title "Velocity profile near wall at x=0.546m"
set xlabel "u, m/s"
set ylabel "y, mm"
set key left top
set yrange [0:20.0]
plot "coles-ix21.dat" using (\$6):((0.24-\$2)*1000) title "Elmer3", \
     "53010801.txt" using (\$7 * 677.4):(\$2*1000) title "Coles" pointtype 4
EOF

gnuplot<<EOF
set term postscript eps 20
set output "coles-x05-temperature.eps"
set title "Temperature profile near wall at x=0.546m"
set xlabel "T, degree K"
set ylabel "y, mm"
set key right top
set yrange [0:20.0]
plot "coles-ix21.dat" using (\$20):((0.24-\$2)*1000) title "Elmer3", \
     "53010801.txt" using (\$8*83.345):(\$2*1000) title "Coles" pointtype 4
EOF

gnuplot<<EOF
set term postscript eps 20
set output "coles-x05-pitot.eps"
set title "Pitot-pressure profile near wall at x=0.546m"
set xlabel "Ppitot, kPa"
set ylabel "y, mm"
set key left top
set yrange [0:20.0]
plot "coles-ix21.dat" using (\$21/1000):((0.24-\$2)*1000) title "Elmer3", \
     "53010801.txt" using (\$3*1.358):(\$2*1000) title "Coles" pointtype 4
EOF

gnuplot<<EOF
set term postscript eps 20
set output "coles-x05-density.eps"
set title "Density profile near wall at x=0.546m"
set xlabel "rho, g/m**3"
set ylabel "y, mm"
set key left top
set yrange [0:20.0]
plot "coles-ix21.dat" using (\$5*1000):((0.24-\$2)*1000) title "Elmer3", \
     "53010801.txt" using (\$9*1000.0*1358.0/(287.0*83.345)/(\$7+0.000001)):(\$2*1000) \
                    title "Coles" pointtype 4
EOF

gnuplot<<EOF
set term postscript eps 20
set output "coles-x05-tke.eps"
set title "Turbulent KE profile near wall at x=0.546m"
set xlabel "tke, m^2/s^2"
set ylabel "y, mm"
set yrange [0:20.0]
plot "coles-ix21.dat" using (\$16):((0.24-\$2)*1000) title "Elmer3"
EOF

gnuplot<<EOF
set term postscript eps 20
set output "coles-x05-omega.eps"
set title "omega profile near wall at x=0.546m"
set xlabel "omega/1.0e6, 1/s"
set ylabel "y, mm"
set yrange [0:20.0]
plot "coles-ix21.dat" using (\$17/1.0e6):((0.24-\$2)*1000) title "Elmer3"
EOF

gnuplot<<EOF
set term postscript eps 20
set output "coles-x05-mu-t.eps"
set title "Turbulence viscosity profile near wall at x=0.546m"
set xlabel "mu-turb, Pa.s"
set ylabel "y, mm"
set yrange [0:20.0]
set xtic 0.0002
plot "coles-ix21.dat" using (\$13):((0.24-\$2)*1000) title "Elmer3"
EOF
