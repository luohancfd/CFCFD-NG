# t4_9945_plot.sh

echo "Recale the pressures and times to get more convenient units."
awk '$1 != "#" {print $1*1.0e3 -259.58, $6/1.0e6}' t4_9945_hx4.dat > pstag_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -259.58, $9/1.0e6}' t4_9945_hx4.dat > hstag_sim.MJkg
awk '$1 != "#" {print $1*1.0e3 -259.58, $8}' t4_9945_hx4.dat > Tstag_sim.degreeK
awk '$1 != "#" {print $1/1.0e3 -6.455, $2/1.0e3}' T4-9945-spa.us-kPa > pstag_exp.MPa
awk '$1 != "#" {print $1*1.0e3 -259.58, $10/1.0e3}' t4_9945_hx5.dat > pitot_sim.kPa
awk '$1 != "#" {print $1*1.0e3 -259.58, $6/1.0e6}' t4_9945_hx1.dat > pshock_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -256.0, $6/1.0e6}' t4_9945_hx0.dat > pcomp_sim.MPa

echo "Pitot pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_9945_ppitot.eps"
set title "T4 Shot 9945: Pitot Pressure"
set xlabel "t, ms"
set ylabel "p, kPa"
set xrange [-1.0:6.0]
plot "pitot_sim.kPa" using 1:2 title "L1d Simulation" with lines linestyle 1
EOF

echo "Nozzle supply pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_9945_pstag.eps"
set title "T4 Shot 9945: Nozzle Supply Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-1.0:6.0]
set yrange [0:40.0]
plot "pstag_sim.MPa" using 1:2 title "L1d Simulation" with lines linestyle 1, \
     "pstag_exp.MPa" using 1:2 title "Measured" with lines linestyle 2
EOF

echo "Nozzle supply temperature"
gnuplot <<EOF
set term postscript eps 20
set output "t4_9945_Tstag.eps"
set title "T4 Shot 9945: Nozzle Supply Temperature"
set xlabel "t, ms"
set ylabel "T, K"
set xrange [-1.0:6.0]
set yrange [0:5000.0]
plot "Tstag_sim.degreeK" using 1:2 title "L1d Simulation" with lines linestyle 1
EOF

echo "Nozzle supply enthalpy"
gnuplot <<EOF
set term postscript eps 20
set output "t4_9945_hstag.eps"
set title "T4 Shot 9945: Nozzle Supply Enthalpy"
set xlabel "t, ms"
set ylabel "H, MJ/kg"
set xrange [-1.0:6.0]
set yrange [0:20.0]
plot "hstag_sim.MJkg" using 1:2 title "L1d Simulation" with lines linestyle 1
EOF

echo "Compression tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_9945_pcomp.eps"
set title "T4 Shot 9945: Compression Tube Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-15.0:10.0]
set yrange [0:40.0]
set key top left
plot "pcomp_sim.MPa" using 1:2 title "L1d Simulation" with lines linestyle 1
EOF

echo "Shock tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_9945_pshock.eps"
set title "T4 Shot 9945: Shock Tube Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-4.0:2.0]
set yrange [0:5.0]
set key bottom right
plot "pshock_sim.MPa" using 1:2 title "L1d Smulation" with lines linestyle 1
EOF
