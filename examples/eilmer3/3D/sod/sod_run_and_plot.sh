# sod_run_and_plot.sh
# Sod's 1-D shock tube exercise as a 3D simulation overkill.
#
e3prep.py --job=sod
time e3shared.exe --job=sod --run
e3post.py --job=sod --output-file=sod_new.dat --slice-list="0:1,:,0,0"

gnuplot<<EOF
set term postscript eps
set output "sod_p.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Pressure, Pa"
set xrange [0.0:1.0]
set yrange [0.0:120.0e3]
plot "sod_new.dat" using 1:9 with points 1 1, \
     "sod_old.dat" using 1:7 with points 1 2
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_rho.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Density, kg/m**3"
set xrange [0.0:1.0]
set yrange [0.0:1.2]
plot "sod_new.dat" using 1:5 with points 1 1, \
     "sod_old.dat" using 1:3 with points 1 2
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_u.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Velocity, m/s"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "sod_new.dat" using 1:6 with points 1 1, \
     "sod_old.dat" using 1:4 with points 1 2
EOF

gnuplot<<EOF
set term postscript eps
set output "sod_T.eps"
set title "One-D Shock Tube at t = 0.6ms"
set xlabel "x, m"
set ylabel "Temperature, K"
set xrange [0.0:1.0]
set yrange [0.0:500.0]
plot "sod_new.dat" using 1:20 with points 1 1, \
     "sod_old.dat" using 1:10 with points 1 2
EOF
