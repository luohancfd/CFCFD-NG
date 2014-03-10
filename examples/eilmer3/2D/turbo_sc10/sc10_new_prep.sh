#! /bin/sh
# sc10_new_prep.sh
e3prep.py --job=sc10_new --do-svg

# Extract the initial solution data and reformat.
e3post.py --job=sc10_new --tindx=0 --vtk-xml

echo At this point, we should be ready to start the simulation.
