#!/bin/bash
# plot_heat_transfer.sh

# Get the heat-flux data around the surface of the sphere.
# These come from blk-1-0 (block 2) and blk-1-1 (block 3)

stages="0 1 2 3 4"
for STAGE in ${stages} 
do
    echo "Stage $STAGE:"
    e3post.py --job=sphere${STAGE} --tindx=9999 --heat-flux-list="2:3,1,-1,:,0" \
        --output-file=sphere_heat_transfer_${STAGE}.dat
done

# Scale current physical simulation to compare with theory and experiment.
awk '$1 != "#" {print $1/0.0066*180.0/3.14159, $2/2.217e6}' \
    sphere_heat_transfer_4.dat > sphere_normalised_heat_transfer.dat

gnuplot <<EOF
set term postscript eps enhanced 20
set output "sphere_norm_heat_transfer.eps"
set style line 1 linetype 1 linewidth 3.0
set xlabel "angle from stagnation point, degrees"
set ylabel "q/q_s"
set logscale y
# set yrange [0.1:2.0]
set yrange [0.1:1.2]
set title "Normalised heat transfer to R6.6mm sphere with Ms=8"
# set key top right
set key bottom left
plot "sphere_normalised_heat_transfer.dat" using 1:2 \
          title "Eilmer3 simulation" with lines ls 1, \
     "kemp_theory.dat" using 1:2 title "Kemp-Rose-Detra theory" \
          with linespoints, \
     "kemp_experiment.dat" using 1:2 title "Kemp-Rose-Detra experiment" \
          with points
EOF

