#! /bin/sh
# 1. Prepare, run and post-process zero blowing case for TINA comparison of heat flux
# for X in mdot_0.000 mdot_0.000-PI
# cd $X
# e3prep.py --job=sphere.py
# mpirun -np 16 e3mpi.exe -f sphere -r
# e3post.py --job=sphere --tindx=9999 --vtk-xml --heat-flux-list="12:15,1,:,:,:"
# cd ..
# done

# 2. Plot
gnuplot << EOF
set term postscript eps enhanced color 'Helvetica,20'
set output 'e3vTINA_heat_flux.eps'
set grid
set size 1.25,1.0
set xlabel "Distance along surface, s (cm)"
set ylabel "Heat transfer rate, q (W/cm**2)"
plot "mdot_0.000/hf_profile.data" u (\$1*100.0):(\$2/1.0e4) w l lw 4 t 'eilmer3: 60 x 60 grid', \
     "mdot_0.000-PI/hf_profile.data" u (\$1*100.0):(\$2/1.0e4) w l lw 4 t 'eilmer3 (PI): 60 x 60 grid', \
     "TINA_solutions/TINA_kaattari_hemisphere_heating.txt" u 1:(\$2*83.0) w l lw 2 lt -1 t 'TINA Solution'
EOF
