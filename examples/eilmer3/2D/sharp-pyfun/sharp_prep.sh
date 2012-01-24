#! /bin/sh
# sharp_prep.sh
# A sharp axisymmetric body as described in Andersons Hypersonics text.

e3prep.py --job=sharp --do-svg

# Extract the initial solution data and reformat so that we can plot the grid.
e3post.py --job=sharp --tindx=0 --vtk-xml

echo At this point, we should be ready to start the simulation.

