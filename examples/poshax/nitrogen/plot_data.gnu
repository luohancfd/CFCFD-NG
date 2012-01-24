# plot_data.gnu

T_inf = 298.0

set term postscript eps enhanced "Palatino-Roman" 20
set output "profile_T.eps"
set xlabel "distance, cm"
set ylabel "temperature, K"
set xrange [0:2]
plot "output-1T.data" u ($1*100.0):2 t 'single temperature' w l 1, \
     "output-2T.data" u ($1*100.0):2 t '2-T: transrotational' w l 2, \
     "output-2T.data" u ($1*100.0):15 t '3-T: vibrational/electronic' w l 3

u_s1 = 706.139
rho_s1 = 0.0136616

u_s2 = 1082.789
rho_s2 = 0.0089093489

set output "profile_u_rho.eps"
set ylabel "velocity, density ratio"
set yrange [0:3]
plot "output-1T.data" u ($1*100.0):($7/u_s1) t '1-T: u ratio' w l 1, \
     "output-1T.data" u ($1*100.0):($4/rho_s1) t '1-T: {/Symbol r} ratio' w l 2, \
     "output-2T.data" u ($1*100.0):($7/u_s2) t '2-T: u ratio' w l 3, \
     "output-2T.data" u ($1*100.0):($4/rho_s2) t '2-T: {/Symbol r} ratio' w l 4

set output "profile_moles.eps"
set ylabel "mole fraction"
set logscale y
set yrange [1.0e-5:1.0]
plot "output-2T.data" u ($1*100.0):9 t '2-temperature model' w l 1, \
     "output-2T.data" u ($1*100.0):10 notitle w l 1, \
     "output-2T.data" u ($1*100.0):11 notitle w l 1, \
     "output-2T.data" u ($1*100.0):12 notitle w l 1, \
     "output-2T.data" u ($1*100.0):13 notitle w l 1, \
     "output-1T.data" u ($1*100.0):9 t '1-temperature model' w l 2, \
     "output-1T.data" u ($1*100.0):10 notitle w l 2, \
     "output-1T.data" u ($1*100.0):11 notitle w l 2, \
     "output-1T.data" u ($1*100.0):12 notitle w l 2, \
     "output-1T.data" u ($1*100.0):13 notitle w l 2

