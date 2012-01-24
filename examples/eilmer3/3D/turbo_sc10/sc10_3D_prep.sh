#! /bin/sh
# sc10_3D_prep.sh
e3prep.py --job=sc10_3D --do-vrml

# Extract the initial solution data and reformat.
e3post.py --job=sc10_3D --tindx=0 --vtk-xml

echo At this point, we should be ready to start the simulation.
