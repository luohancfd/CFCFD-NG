#! /bin/sh
# run_extended.sh
# exercise the Navier-Stokes solver for the cone20-parametric test case.
# It is assumed that the path is set correctly.

# Prepare the simulation input files (parameter, grid and initial flow data).
# The SVG file provides us with a graphical check on the geometry
e3prep.py --job=conepe --do-svg
if [ "$?" -ne "0" ] ; then
    echo "e3prep.py ended abnormally."
    exit
fi

# Integrate the solution in time, 
# recording the axial force on the cone surface.
time e3shared.exe -f conepe --run --verbose
if [ "$?" -ne "0" ] ; then
    echo "e3shared.exe ended abnormally."
    exit
fi

# Extract the solution data and reformat.
# If no time is specified, the final solution found is output.
e3post.py --job=conepe --vtk-xml --add-mach

echo "At this point, we should have a solution that can be viewed with Paraview."
