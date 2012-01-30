#!/bin/bash
# post_simulation.sh

# Create a VTK plot file of the steady full flow field.
e3post.py --job=cyl --tindx=9999 --vtk-xml

# Pull out the cylinder surfaces.
e3post.py --job=cyl --tindx=9999 --output-file=cylinder \
    --surface-list="0,east;1,east;3,bottom"

# Now pull out some block surfaces that show cross-sections of the flow field.
e3post.py --job=cyl --tindx=9999 --output-file=interior \
    --surface-list="0,bottom;1,bottom;0,north;1,north;2,north;3,north;0,south;1,south;2,south;3,south;3,east"

# Stagnation-line flow data
e3post.py --job=cyl --tindx=9999 --slice-list="0,:,0,0" \
    --output-file=stagnation-line.data

# Plot temperature and N2 density profiles along the stagnation-line
# NOTE: thermal equilibrium solution needs to be present
gnuplot <<EOF
set term postscript eps enhanced "Helvetica" 20
set output "temperature_profiles.eps"
set size 1.0,1.0
set ylabel "Temperature (K)"
set xlabel "Distance from stagnation point, x (mm)"
set grid
set key top right
plot '../thermal-eq/stagnation-line.data' u (\$1*1000+7.5):25 w l lt 1 lw 3 t "Thermal eq.: T", \
     'stagnation-line.data' u (\$1*1000+7.5):25 w lp lt 2 lw 2 pt 4 ps 0.7 t "Thermal noneq.: T_{tr}", \
     'stagnation-line.data' u (\$1*1000+7.5):27 w lp lt 3 lw 2 pt 5 ps 0.7 t "Thermal noneq.: T_{ve}"

set output "N2_profiles.eps"
set key top right
set ylabel "N_2 density (kg/m^3)"
plot '../thermal-eq/stagnation-line.data' u (\$1*1000+7.5):(\$5*\$18) w l lt 1 lw 3 t "Thermal eq.", \
     'stagnation-line.data' u (\$1*1000+7.5):(\$5*\$18) w lp lt 2 lw 2 pt 4 ps 0.7 t "Thermal noneq."
EOF
epstopdf temperature_profiles.eps
epstopdf N2_profiles.eps

