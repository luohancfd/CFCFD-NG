#!/bin/bash
# post_simulation.sh

# Extract the stagnation line data from the steady flow field.
e3post.py --job=n90 --output-file=n90_100_iy1.data --tindx=5 \
    --slice-list="0,:,1,0"
gnuplot plot_comparison.gnu

# Create a VTK plot file of the steady flow field.
e3post.py --job=n90 --tindx=5 --vtk-xml

