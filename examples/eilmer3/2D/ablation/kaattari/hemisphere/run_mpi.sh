#! /bin/sh
# 1. Prepare, run and post-process all cases (starting with mdot=0 to enable others to use this solution)
for D in mdot_0.000 mdot_0.004 mdot_0.010 mdot_0.013
do
cd $D
e3prep.py --job=sphere.py
mpirun -np 16 e3mpi.exe -f sphere -r
e3post.py --job=sphere --tindx=9999 --vtk-xml --heat-flux-list="12:15,1,:,:,:"
cd ..
done

# 2. Plot
gnuplot << EOF
set term postscript eps enhanced color 'Helvetica,20'
set output 'heating_profiles.eps'
set grid
set 
set size 1.25,1.0
set xlabel "Distance along surface, s (cm)"
set ylabel "Heat transfer rate, q/q_0"
plot "mdot_0.000/hf_profile.data" u ($1*100.0):($2/8.307181e+05) w l lw 4 t 'eilmer3: mdot* = 0.0', \
     "mdot_0.004/hf_profile.data" u ($1*100.0):($2/8.307181e+05) w l lw 4 t 'eilmer3: mdot* = 0.004', \
     "mdot_0.010/hf_profile.data" u ($1*100.0):($2/8.307181e+05) w l lw 4 t 'eilmer3: mdot* = 0.010', \
     "mdot_0.013/hf_profile.data" u ($1*100.0):($2/8.307181e+05) w l lw 4 t 'eilmer3: mdot* = 0.013', \
     "kaattari_data/mdot_0.004.txt" u 1:2 w p ps 2 pt 4 lt 2 t 'Kaattari: mdot* = 0.004', \
     "kaattari_data/mdot_0.010.txt" u 1:2 w p ps 2 pt 5 lt 3 t 'Kaattari: mdot* = 0.010', \
     "kaattari_data/mdot_0.013.txt" u 1:2 w p ps 2 pt 6 lt 4 t 'Kaattari: mdot* = 0.013'
EOF
