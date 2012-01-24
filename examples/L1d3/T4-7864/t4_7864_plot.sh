# t4_7864_plot.sh

echo "Recale the pressures and times to get more convenient units."
awk '$1 != "#" {print $1*1.0e3 -149.77, $6/1.0e6}' t4_7864_hx4.dat > pstag_sim.MPa
awk '$1 != "#" {print $1/1.0e3 -6.455, $2/1.0e3}' T4-7864-spa.us-kPa > pstag_exp.MPa
awk '$1 != "#" {print $1*1.0e3 -149.77, $10/1.0e3}' t4_7864_hx5.dat > pitot_sim.kPa
awk '$1 != "#" {print $1*1.0e3 -149.77, $6/1.0e6}' t4_7864_hx1.dat > pshock_sim.MPa
awk '$1 != "#" {print $1*1.0e3 -144.17, $6/1.0e6}' t4_7864_hx0.dat > pcomp_sim.MPa

echo "Pitot pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_7864_ppitot.eps"
set title "T4 Shot 7864: Pitot Pressure"
set xlabel "t, ms"
set ylabel "p, kPa"
set xrange [-1.0:6.0]
plot "pitot_sim.kPa" using 1:2 title "L1d Simulation" with lines linestyle 1
EOF

echo "Nozzle supply pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_7864_pstag.eps"
set title "T4 Shot 7864: Nozzle Supply Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-1.0:6.0]
set yrange [0:100.0]
plot "pstag_sim.MPa" using 1:2 title "L1d Simulation" with lines linestyle 1, \
     "pstag_exp.MPa" using 1:2 title "Measured" with lines linestyle 2
EOF

echo "Compression tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_7864_pcomp.eps"
set title "T4 Shot 7864: Compression Tube Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-15.0:10.0]
set yrange [0:100.0]
set key top left
plot "pcomp_sim.MPa" using 1:2 title "L1d Simulation" with lines linestyle 1
EOF

echo "Shock tube pressure"
gnuplot <<EOF
set term postscript eps 20
set output "t4_7864_pshock.eps"
set title "T4 Shot 7864: Shock Tube Pressure"
set xlabel "t, ms"
set ylabel "p, MPa"
set xrange [-4.0:2.0]
set yrange [0:20.0]
set key bottom right
plot "pshock_sim.MPa" using 1:2 title "L1d Smulation" with lines linestyle 1
EOF
