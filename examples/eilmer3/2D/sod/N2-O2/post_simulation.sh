#!/bin/bash

# Extract the profile along the shock tube.
# slice-list=block-range,i-range,j-range,k-range
#            all blocks - :
#            all i's - :
#            constant j - 0
#            constant k - 0 (not relevant in 2D anyway)
e3post.py --job=sod --slice-list=":,:,0,0" --output-file=profile.data

# Plot with Gnuplot
gnuplot plot_profile.gplot


