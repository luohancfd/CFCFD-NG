#! /bin/sh
# post.sh

JOB=lam_flat_plate
TINDX=9999

# Extracts Paraview or Visit files for flowfield visualisation
e3post.py --job=$JOB --vtk-xml --tindx=$TINDX

# Extracts slices at x = 1.0 m. 
# This will allow us to examine the boundary layer profile.
e3post.py --job=$JOB --output-file=profile-at-x1m.dat --tindx=$TINDX \
    --slice-at-point="12,jk,1.0,0.1,0.0;13,jk,1.0,0.1,0.0;14,jk,1.0,0.1,0.0;15,jk,1.0,0.1,0.0" \
    --add-pitot-p --add-mach

awk -f integral-thicknesses.awk profile-at-x1m.dat > thicknesses.txt

gnuplot<<EOF
set term postscript eps enhanced 20
set output "u-velocity-comparison.eps"
set title "Velocity profile near wall at x = 1.0 m"
set xlabel "u/u_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./profile-at-x1m.dat" using ((\$6)/1390.0):((0.44-\$2)*1000) \
     title "Eilmer3" with points pt 1, \
     "boundary-layer-profile-clbl.data" using (\$3):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF

gnuplot<<EOF
set term postscript eps enhanced 20
set output "density-comparison.eps"
set title "Density profile near wall at x = 1.0 m"
set xlabel "rho/rho_e"
set ylabel "y, mm"
set key left top
set xrange [0:1.2]
set yrange [0:20]
plot "./profile-at-x1m.dat" using ((\$5)/0.011835):((0.44-\$2)*1000) \
     title "Eilmer3" with points pt 1, \
     "boundary-layer-profile-clbl.data" using (\$5):((\$2)*1000) \
     title "CLBL" with lines lt 1
EOF

# Extracts a slice of the nearest cells (to the wall) along the plate. 
# This will allow us to examine viscous properties, like skin-friction 
# and y+ values.
e3post.py --job=$JOB --output-file=first-cell-off-wall.dat --tindx=$TINDX \
    --slice-list="3,:,-1,0;7,:,-1,0;11,:,-1,0;13,:,-1,0" \
    --add-pitot-p --add-mach

python compute_shear_stress.py
python reference_temperature_method.py

gnuplot<<EOF
set term postscript eps enhanced 20
set output "cf-comparison.eps"
set title "Skin friction coefficient along plate"
set xlabel "x, m"
set ylabel "c_f"
set key left top
set yrange [0:0.004]
set key top right
plot "cf-eilmer3.data" using (\$1):(\$5) title "Eilmer3" with points pt 1, \
     "cf-clbl.data" using (((\$1)-1)/200.0):(\$2) title "CLBL" with lines lt 1, \
     "cf-ref-temp.data" using (\$1):(\$2) title "RefTemp" with lines lt 2
EOF
