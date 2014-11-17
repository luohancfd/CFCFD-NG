#! /bin/sh
# Kelvin_Helmholtz_run_mpi.sh
# exercise the Navier-Stokes solver with added MHD capability for the MHD-test test case.
# It is assumed that the path is set correctly.

### do the non-clean case first:

cd ./un-clean

# Prepare the simulation input files (parameter, grid and initial flow data).
e3prep.py --job=Kelvin-Helmholtz_no_clean --do-svg
if [ "$?" -ne "0" ] ; then
    echo "e3prep.py ended abnormally."
    exit
fi

# Integrate the solution in time 
mpirun -np 4 e3mpi.exe -f Kelvin-Helmholtz_no_clean --run --verbose
if [ "$?" -ne "0" ] ; then
    echo "e3shared.exe ended abnormally."
    exit
fi

# Extract the solution data and reformat.
e3post.py --job=Kelvin-Helmholtz_no_clean --vtk-xml --tindx=all

cd ..

### now do the clean case:

cd ./clean

# Prepare the simulation input files (parameter, grid and initial flow data).
e3prep.py --job=Kelvin-Helmholtz_clean --do-svg
if [ "$?" -ne "0" ] ; then
    echo "e3prep.py ended abnormally."
    exit
fi

# Integrate the solution in time 
mpirun -np 4 e3mpi.exe -f Kelvin-Helmholtz_clean --run --verbose
if [ "$?" -ne "0" ] ; then
    echo "e3shared.exe ended abnormally."
    exit
fi

# Extract the solution data and reformat.
e3post.py --job=Kelvin-Helmholtz_clean --vtk-xml --tindx=all

cd ..

### Post-processing

echo "At this point, we should have a solution that can be viewed with Paraview."

# Plot the divergence of magnetic field over time for the two cases just run
python divB-over-time.py
echo "Should now have a plot of divergence vs time"

