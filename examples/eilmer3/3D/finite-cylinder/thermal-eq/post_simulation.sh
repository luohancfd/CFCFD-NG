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
