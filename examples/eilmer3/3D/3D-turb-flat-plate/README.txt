# README.txt for 3D-turb-flat-plate (turb_flat_plate.py)
# Test Case 1 - Three dimensional flat plate utilising the k-omega turbulence model
# Author: Sam Stennett, Nov 2014

This is one of the test cases for Samuel Stennett's undergraduate thesis (2014).
It is an extension of the 2D example case turb_flat_plate.py

The directory contains the images indicating correlation along the width, 
a grid convergence study with an 80% reduced mesh, as well as a comparison 
with the original 2D case.

There are sufficient files in the directory to replicate the results (apart from
the grid convergence study). Some files have been condensed from original versions
in order to reduce the number of files in the directory. They are expected to work,
and any issues that arise should be easy to solve.

To run the test case, enter the terminal and execute the following commands:

	./prep.sh
	./run.sh
	./post.sh

Note: post.sh contains the execute scripts for plot.sh and plot2Dv3D.sh, so these 
scripts are not required to be executed separately.

The run.sh script is written for execution on the HPCU cluster `Barrine', or any
system utilising PBS queue software and the MPI version of Eilmer. The expected 
runtime of the case is ~1 week, utilising 16 CPUs. Increasing the number of blocks 
and CPUs would speed up the simulation. If no MPI version is installed, or an 
insufficient number of CPUs, the case can instead be run using the command:

	e3shared.exe --job=turb_flat_plate --run

The runtime of the case executed in this manner will be very long (months?) on a 
typical PC or laptop. It is highly recommended that `Barrine' access is granted in
order to run the 3D simulations, due to the large number of cells and use of the 
k-omega turbulence model.

----------------------------------------------------------------------------

