#! /bin/sh
# 1. Prepare, run and post-process all cases (starting with mdot=0 to enable others to use this solution)
for D in mdot_0.000 mdot_0.002 mdot_0.009
do
cd $D
e3prep.py --job=sphere.py
e3shared.exe -f sphere -r
e3post.py --job=sphere --tindx=9999 --vtk-xml --heat-flux-list="6,1,:,:,:;7,1,:,:,:;8,1,:,:,:;11,1,:,:,:"
cd ..
done

# 2. Plot
gnuplot << EOF
set term postscript eps enhanced color 'Helvetica,20'
set output 'heating_profiles.eps'
set size 1.25,1.0
set grid
set xlabel "Distance along surface, s (cm)"
set ylabel "Heat transfer, q (W/cm**2)"
plot "mdot_0.000/hf_profile.data" u ($1*100):($2/8.307181e+05) w lp lw 4 t 'mdot_wall/mdot_inf = 0.000', \
     "mdot_0.002/hf_profile.data" u ($1*100):($2/8.307181e+05) w lp lw 4 t 'mdot_wall/mdot_inf = 0.002', \
     "mdot_0.009/hf_profile.data" u ($1*100):($2/8.307181e+05) w lp lw 4 t 'mdot_wall/mdot_inf = 0.009', \
     "kaattari_data/mdot_0.002.txt" u 1:2 w p ps 2 pt 4 lt 2 t 'Kaattari: mdot* = 0.002', \
     "kaattari_data/mdot_0.009.txt" u 1:2 w p ps 2 pt 5 lt 3 t 'Kaattari: mdot* = 0.009'
EOF
