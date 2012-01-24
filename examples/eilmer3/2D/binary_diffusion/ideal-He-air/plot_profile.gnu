

#set label 1 "f_{He}" at -15,0.8
#set label 2 "f_{air}" at 15,0.8
#set arrow 1 from -13.8,0.82 to -10,0.91
#set arrow 2 from 14,0.82 to 10,0.91
set term postscript eps enhanced 20
set output "mass_f_profile.eps"
set title "Binary diffusion of helium and air\np= 1 atm, T= 300.0K, t=5.0e-7s"
set xlabel "r, {/Symbol m}m"
set ylabel "mass fraction"
set key 16,0.6
set yrange [-0.1:1.1]
plot "numerical-profile-Ficks.data" u ($1*1e6):18 t "Fick's first law" w p 1, \
     "numerical-profile-Ficks.data" u ($1*1e6):19 notitle w p 1
