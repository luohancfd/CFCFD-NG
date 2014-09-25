#! /bin/sh
# cone20_run_mpi.sh
# exercise the Navier-Stokes solver for the Cone20 test case.
# It is assumed that the path is set correctly.

# Prepare the simulation input files (parameter, grid and initial flow data).
# The SVG file provides us with a graphical check on the geometry
e3prep.py --job=cone20 --do-svg

# Integrate the solution in time using the parallel processing code. 
mpirun -np 2 e3mpi.exe -f cone20 --run

# Extract the solution data and reformat.
# If no time is specified, the final solution found is output.
e3post.py --job=cone20 --vtk-xml

echo "At this point, we should have a solution that can be viewed with Paraview."

