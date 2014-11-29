# README.txt for mallinson_cylinder
# Author: Sam Stennett, Nov 2014

This is one of the test cases for the k-omega validation report.
The case was written by Wilson Y. K. Chan in 2008 and edited for the repository 
by Samuel Stennett in November 2014.

The directory contains the images for case including experimental setup, an explanation
of the turbulence-zone option, y+ information on sensitvity and a typical distribution,
as well as heat flux and pressure data comparisons of differing y+ and aspect ratio, 
boundary layer profiles at 800mm and the final solutions vs experimental data. 

There are sufficient files in the directory to replicate the results. 
Some files have been condensed from original versionsin order to reduce the number 
of files in the directory.


################################ CASE PREPARATION ################################


To run the test case, first:

	SELECT A DESIRED RESOLUTION! 

..by specifying the 'casenum' variable (from 0 to 8) inside the mallinson_cylinder.py 
file. Once selected, enter the terminal and execute:

	e3prep.py --job=mallinson_cylinder --verbosity=1
	

################################# CASE EXECUTION #################################


Depending on the the case selected, the following run command should be executed:

 - FOR CASE 0 (100x33cells-yplus3.0-ar215):
	
	./run4.sh

 - FOR CASE 1 & 2 (200x65cells cases):

	./run8.sh

 - FOR CASES 3-8 (400x130cells, 600x130cells, 520x169cells cases):

	./run16.sh

It should be noted that the walltimes that have been set for the scripts is the
maximum for cases using the same number of blocks/CPUs. This number can be
reduced depending on the number of cells/clustering, as required.

The run.sh script is written for execution on the HPC cluster `Barrine', or any
system utilising PBS queue software and the MPI version of Eilmer. The expected 
runtime of the case depends on the grid resolution and clustering. 
If no MPI version is installed, or an insufficient number of CPUs, the case can 
instead be run using the command:

	e3shared.exe --job=mallinson_cylinder --run

The runtime of the case executed in this manner will be very long (months?) on a 
typical PC or laptop. It is highly recommended that `Barrine' access is granted in
order to run the simulation.


############################## CASE POST-PROCESSING ##############################


The post_EditPerCase.sh script should be edited as required in order to extract the
correct data from the cells closest to the x and y axis (the walls). The block info
and cell info is all that needs to be altered (unless a different time index is 
desired). Following any alterations, the script can be executed using the command:

	./post_EditPerCase.sh

----------------------------------------------------------------------------------
