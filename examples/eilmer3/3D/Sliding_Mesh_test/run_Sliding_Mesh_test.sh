#! /bin/sh
# run_Sliding_Mesh_test.sh
# Generate example data for radial inflw turbine Rotor meshing example
# It is assumed that the path is set correctly.

# Create the grid.
e3prep.py --job=Sliding_Mesh_test
if [ "$?" -ne "0" ] ; then
    echo "e3prep.py ended abnormally."
    exit
fi


# run the code
e3shared.exe --job=Sliding_Mesh_test --run

# Reformat to vtk.
e3post.py --job=Sliding_Mesh_test --vtk-xml --tindx=all --add-mach

# Open paraview
paraview plot/Sliding_Mesh_test.pvd


