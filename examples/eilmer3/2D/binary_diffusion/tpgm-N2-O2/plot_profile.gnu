

set label 1 "f_{N_2}" at -15,0.8
set label 2 "f_{O_2}" at 15,0.8
set arrow 1 from -13.8,0.82 to -10,0.91
set arrow 2 from 14,0.82 to 10,0.91
set term postscript eps enhanced 20
set output "mass_f_profile.eps"
set title "Binary diffusion of nitrogen and oxygen\np= 1 atm, T= 273.2K, t=1.0e-6s"
set xlabel "r, {/Symbol m}m"
set ylabel "mass fraction"
set key 19,0.6
set yrange [-0.1:1.1]
plot "numerical-profile-Ficks.data" u ($1*1e6):18 t "Fick's first law" w p 1, \
     "numerical-profile-Ficks.data" u ($1*1e6):19 notitle w p 1, \
     "exact-profile.data" u ($1*1e6):3 t 'exact solution' w l 2, \
     "exact-profile.data" u ($1*1e6):4 notitle w l 2
