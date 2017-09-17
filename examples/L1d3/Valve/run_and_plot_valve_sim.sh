#! /bin/sh
# run_and_plot_valve_sim.sh


echo 
echo Beginning simulation...
echo Console output is being caught in file valve.log.
time l1d.exe -f valve > valve.log
echo Finished simulation.
echo 


sptime.exe -f valve -tstop 10.0e-2 -p -log
mb_cont.exe -fi valve_log_p.gen -fo valve_log_p.ps -ps \
            -var 2 -edge -notrueshape -xrange 0.0 10.0 1.0 \
            -yrange 0.0 2.0e-2 1.0e-2 -colour


# (3) Flow history at a specified location can be extracted from
#     the history file piston.hx and plotted with GNU-Plot.
#     Since one history location was specified in the parameter file,
#     the only available history location has an index of 0 and its
#     data will be written to the file piston_hx0.dat.

l_hist.exe -f valve -tstart 0.0 -tstop 20.0e-2 -xloc 0
l_hist.exe -f valve -tstart 0.0 -tstop 20.0e-2 -xloc 1
l_hist.exe -f valve -tstart 0.0 -tstop 20.0e-2 -xloc 2
l_hist.exe -f valve -tstart 0.0 -tstop 20.0e-2 -xloc 3


awk '$1 != "#" {print $1 * 1000.0,$3*$4}' valve_hx0.dat > mdot0.dat
awk '$1 != "#" {print $1 * 1000.0,$3 * $4}' valve_hx1.dat > mdot1.dat
awk '$1 != "#" {print $1 * 1000.0,$3 * $4}' valve_hx2.dat > mdot2.dat
awk '$1 != "#" {print $1 * 1000.0,$3 * $4}' valve_hx3.dat > mdot3.dat


gnuplot <<EOF
set term postscript eps 20
set output "pressure_valve_0.eps"
set title "Test Valve"
set xlabel "t, s"
set ylabel "p, Pa"
set xrange [0:0.02]
set key top right
plot "valve_hx0.dat" using 1:6 title "End wall" with lines linestyle 1
EOF

gnuplot <<EOF
set term postscript eps 20
set output "pressure_valve_1.eps"
set title "Test Valve"
set xlabel "t, s"
set ylabel "p, Pa"
set xrange [0:0.02]
set xtic 0.010
set key top right
plot "valve_hx1.dat" using 1:6 title "End wall" with lines linestyle 1
EOF

gnuplot <<EOF
set term postscript eps 20
set output "pressure_valve_2.eps"
set title "Test Valve"
set xlabel "t, s"
set ylabel "p, Pa"
set xrange [0:0.02]
set xtic 0.010
set key top right
plot "valve_hx2.dat" using 1:6 title "End wall" with lines linestyle 1
EOF

gnuplot <<EOF
set term postscript eps 20
set output "pressure_valve_3.eps"
set title "Test Valve"
set xlabel "t, s"
set ylabel "p, Pa"
set xrange [0:0.02]
set xtic 0.010
set key top right
plot "valve_hx3.dat" using 1:6 title "End wall" with lines linestyle 1
EOF

gnuplot <<EOF
set term postscript eps 20
set output "pressure_valve_4.eps"
set title "Test Valve"
set xlabel "t, s"
set ylabel "p, Pa"
set xrange [0:0.02]
set xtic 0.010
set key top right
plot "valve_hx4.dat" using 1:6 title "End wall" with lines linestyle 1
EOF
