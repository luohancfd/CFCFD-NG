#! /bin/sh
# post.sh

JOB=lam_flat_plate_march
TINDX=last

# Extracts Paraview or Visit files for flowfield visualisation
e3post.py --job=$JOB --vtk-xml --tindx=$TINDX

# Extracts slices at x = 1.0 m. 
# This will allow us to examine the boundary layer profile.
e3post.py --job=$JOB --output-file=profile-at-x1m-march.dat --tindx=$TINDX \
    --slice-at-point="84,jk,1.0,0.1,0.0;85,jk,1.0,0.1,0.0;86,jk,1.0,0.1,0.0;87,jk,1.0,0.1,0.0" \
    --add-pitot-p --add-mach

awk -f integral-thicknesses.awk profile-at-x1m-march.dat > thicknesses-march.txt

gnuplot<<EOF
set term postscript eps enhanced 20
set output "u-velocity-comparison-march.eps"
set title "Velocity profile near wall at x = 1.0 m"
set xlabel "u/u_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./profile-at-x1m-march.dat" using ((\$6)/1390.0):((0.44-\$2)*1000) \
     title "e3march.py" with points pt 1, \
     "boundary-layer-profile-clbl.data" using (\$3):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20
set output "density-comparison-march.eps"
set title "Density profile near wall at x = 1.0 m"
set xlabel "rho/rho_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./profile-at-x1m-march.dat" using ((\$5)/0.011835):((0.44-\$2)*1000) \
     title "e3march.py" with points pt 1, \
     "boundary-layer-profile-clbl.data" using (\$5):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF

# Extracts a slice of the nearest cells (to the wall) along the plate. 
# This will allow us to examine viscous properties, like skin-friction 
# and y+ values.
e3post.py --job=$JOB --output-file=first-cell-off-wall-march.dat --tindx=$TINDX \
    --slice-list="3,:,-1,0;7,:,-1,0;11,:,-1,0;15,:,-1,0;19,:,-1,0;23,:,-1,0;27,:,-1,0;31,:,-1,0;35,:,-1,0;39,:,-1,0;43,:,-1,0;47,:,-1,0;51,:,-1,0;55,:,-1,0;59,:,-1,0;63,:,-1,0;67,:,-1,0;71,:,-1,0;75,:,-1,0;79,:,-1,0;83,:,-1,0;87,:,-1,0" \
    --add-pitot-p --add-mach

python compute_shear_stress_march.py
python reference_temperature_method.py

gnuplot<<EOF
set term postscript eps enhanced 20
set output "cf-comparison-march.eps"
set title "Skin friction coefficient along plate"
set xlabel "x, m"
set ylabel "c_f"
set key left top
set yrange [0:0.004]
set key top right
plot "cf-eilmer3-march.data" using (\$1):(\$5) title "e3march.py" with points pt 1, \
     "cf-clbl.data" using (((\$1)-1)/200.0):(\$2) title "CLBL" with lines lt 1, \
     "cf-ref-temp.data" using (\$1):(\$2) title "RefTemp" with lines lt 2
EOF
