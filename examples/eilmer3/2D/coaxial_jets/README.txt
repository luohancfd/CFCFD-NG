# README.txt for mallinson_cylinder
# Author: Sam Stennett, Nov 2014

This is one of the test cases for the k-omega validation report.
The case was written by Wilson Y. K. Chan in 2008 and edited for the repository 
by Samuel Stennett in November 2014.

##################################################################################
IMPORTANT: It should be noted that in porting the files over from older versions
of the solver, the k-omega turbulence model has broken for this case. The current
version runs without a turbulence model, instead producing a laminar solution. It
is recommended that the case is investigated in the future in order to determine
what the cause of the breakage is, and how it can be fixed.
##################################################################################

The directory contains the images for case including experimental facility, 
the grid utilised as well as experimental and numerical Schlieren images for
validation purposes.

There are sufficient files in the directory to replicate the results. 
Some files have been condensed from original versions in order to reduce the number 
of files in the directory.

To run the test case, enter the terminal and execute:

	e3prep.py --job=mallinson_cylinder --verbosity=1

Upon completion of the prep, the following command should be executed:
	
	./run.sh

It should be noted that the walltime that has been set for the case is a ballpark
figure, and will probably require extension. 

The run.sh script is written for execution on the HPC cluster `Barrine', or any
system utilising PBS queue software and the MPI version of Eilmer. The expected 
runtime of the case depends on the grid resolution and clustering. 

If no MPI version is installed, or an insufficient number of CPUs, the case can 
instead be run using the command:

	e3shared.exe --job=mallinson_cylinder --run

The runtime of the case executed in this manner will be very long (months?) on a 
typical PC or laptop. It is highly recommended that `Barrine' access is granted in
order to run the simulation.

The post.sh script should be edited as required in order to extract data at the 
desired time (at a specified time index). Following any alterations, the script 
can be executed using the command:

	./post.sh

----------------------------------------------------------------------------------
