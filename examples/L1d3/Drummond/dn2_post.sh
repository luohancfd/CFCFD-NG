#! /bin/sh
# dn2_post.sh
#
# At this point, the simulation is complete and the post-processing
# programs can be used to reformat the data for plotting.
#
# For example:
# (1) A space-time plot can be generated via a 
#     contour plotting package.
#     First, extract the pressure (on a logarithmic scale) and 
#     save it in a file (dn2_log_p.gen) suitable for 
#     the Generic contour  plotter.
#     Then contour variable 2 (0 is x, 1 is time) to get 
#     the X-T diagram.

sptime.exe -f dn2 -tstop 8.0e-3 -p -log
mb_cont.exe -fi dn2_log_p.gen -fo dn2_log_p.ps -ps \
            -var 2 -edge -notrueshape -xrange -4.0 0.5 1.0 \
            -yrange 0.0 8.0e-3 2.0e-3 -levels 3.5 6.5 0.1

# (2) Flow history at specified locations can be extracted from
#     the history file dn2.hx and plotted with GNU-Plot.

l_hist.exe -f dn2 -tstart 1.0e-3 -tstop 8.0e-3 -xloc 0
l_hist.exe -f dn2 -tstart 2.0e-3 -tstop 8.0e-3 -xloc 1
l_hist.exe -f dn2 -tstart 2.0e-3 -tstop 8.0e-3 -xloc 2
l_hist.exe -f dn2 -tstart 2.0e-3 -tstop 8.0e-3 -xloc 3
l_hist.exe -f dn2 -tstart 2.0e-3 -tstop 8.0e-3 -xloc 4

# Recale the pressures and heat-transfer to get more convenient units.
awk '$1 != "#" {print $1 * 1000.0, $6 / 1.0e6}' dn2_hx1.dat > pstag.Mpa
awk '$1 != "#" {print $1 * 1000.0, $10 / 1.0e3}' dn2_hx4.dat > pitot.kpa
awk '$1 != "#" {print $1 * 1000.0, $12 / -1.0e6}' dn2_hx2.dat > heatflux.Mwatt

gnuplot <<EOF
# was file ppitot.gnu
set term postscript eps 20
set output "dn2_ppitot.eps"
set title "Drummond Tunnel with N2 Driver, Pitot Pressure"
set xlabel "t, ms"
set ylabel "p, kPa"
set xrange [3.0:8.0]
set xtic 1.0
set key top left
plot "pitot.kpa" using 1:2 title "L1d Simulation" with lines linestyle 1
EOF

gnuplot <<EOF 
# was file heat_trans.gnu
set term postscript eps 20
set output "dn2_heatflux.eps"
set title "Drummond Tunnel with N2 Driver, Heat Transfer"
set xlabel "t, ms"
set ylabel "Q, MW/m**2"
set xrange [3.0:8.0]
set xtic 1.0
set key top left
plot "heatflux.Mwatt" using 1:2 title "Simulation" with lines linestyle 1
EOF

gnuplot <<EOF
# was file pstag.gnu
set term postscript eps 20
set output "dn2_pstag.eps"
set title "Drummond Tunnel with N2 Driver, Nozzle Supply Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [3.0:8.0]
set xtic 1.0
set yrange [0:2.0]
set key top left
plot "pstag.Mpa" using 1:2 title "Simulation" with lines linestyle 1
EOF

# Now set up a couple of files that describe the tube cross-section.
awk '$1 != "#" {print $1, $2 / 2.0}' dn2.dump > upper.profile
awk '$1 != "#" {print $1, -$2 / 2.0}' dn2.dump > lower.profile

gnuplot <<EOF
# was file d_tube.gnu
set term postscript eps 20
set output "dn2_tube.eps"
set title "Drummond Tunnel Facility Profile"
set xlabel "x, m"
set ylabel "y, m"
set xrange [-4:0.5]
set yrange [-0.1:0.1]
plot "upper.profile" using 1:2 title "" with lines linestyle 1, \
     "lower.profile" using 1:2 title "" with lines linestyle 1
EOF

gnuplot <<EOF
# was file d_noz.gnu
set term postscript eps 20
set output "dn2_noz.eps"
set title "Drummond Tunnel Nozzle Profile"
set xlabel "x, m"
set ylabel "y, m"
set xrange [0:0.3]
set yrange [-0.1:0.1]
plot "upper.profile" using 1:2 title "" with lines linestyle 1, \
     "lower.profile" using 1:2 title "" with lines linestyle 1
EOF