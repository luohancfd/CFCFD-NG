#! /bin/sh
# kiock_prep.sh
e3prep.py --job=kiock --do-svg

# Extract the initial solution data and reformat.
e3post.py --job=kiock --tindx=0 --vtk-xml

echo At this point, we should be ready to start the simulation.
