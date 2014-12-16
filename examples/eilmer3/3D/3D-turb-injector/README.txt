# README.txt for 3D-turb-injector (inject.py)
# Test Case 3 - Three dimensional injector utilising the k-omega turbulence model
# Author: Sam Stennett, Dec 2014

This is one of the test cases for Samuel Stennett's undergraduate thesis (2014).

The directory contains the images indicating the basic flow schematic, and the 
block arrangement and mesh used in the 3D simulation. Under the results-figures
folder, graphs of pressure coefficient along the flat plate surface (on the 
symmetry plane) are available, Mach and Cp figures from Paraview, as well as 
the results from Viti and Wallis for use in the validation.

There are sufficient files in the directory to produce a 2D `run-up', interpolate
and extrapolate a 3D inflow profile, and simulate the complete 3D flow domain. 
The case can also be run as is.

################################
#### RUNNING THE CASE AS IS ####
################################

To run the test case, enter the terminal and execute the following commands:

	e3prep.py --job=inject --verbosity=1

The verbosity flag allows for feedback during the prepping process, as the large
number of cells and blocks requires some time to prepare.
Next, the .config file must be edited to input the appropriate inflow profiles. 
If the .py script has NOT been altered, this can be done by executing the following:

	cp inject.config.backup inject.config

This will overwrite the current .config file with the version containing the new
boundary conditions. This can also be executed manually by editing the .config file
and replacing the REPLACE.dat entries with the appropriate .dat files. Currently, 
the geometry has 16 inflow boundaries (BLOCKS 0-15, WEST boundary), and can be set 
in ascending order repetitively:

	for B0: 	3Dinfilesmalllower1.dat
	for B1: 	3Dinfilesmalllower2.dat
	for B2: 	3Dinfilesmalllower3.dat
	for B3: 	3Dinfilesmalllower4.dat
	for B4: 	3Dinfilesmallupper1.dat
	for B5: 	3Dinfilesmallupper2.dat
	for B6: 	3Dinfilesmallupper3.dat
	for B7: 	3Dinfilesmallupper4.dat
	for B8: 	3Dinfilelargelower1.dat
	...		...
	for B15: 	3Dinfilelargeupper4.dat

This is done as the small blocks are defined vertically (0 to 7), and then step sideways
to the large blocks and repeat.
Finally, the case can be run by executing:

	./run_barrine.sh

The run_barrine.sh script is written for execution on the HPCU cluster `Barrine', 
or any system utilising PBS queue software and the MPI version of Eilmer. The expected 
runtime of the case is over a month, utilising 96 CPUs. 
If no MPI version is installed, or an insufficient number of CPUs, the case can 
instead be run using the command:

	e3shared.exe --job=inject --run

The runtime of the case executed in this manner will be very long (years?) on a 
typical PC or laptop. It is highly recommended that `Barrine' access is granted in
order to run the 3D simulations, due to the large number of cells and use of the 
k-omega turbulence model.

Additionally, sufficient post-processing scripts are available in the post-files 
folder to produce plots of pressure coefficient along the flat plate centreline, 
and the associated Paraview files. The files in the post-files folder should be 
copied into the 3D-turb-injector directory before execution - this was done to 
clean up the directory. 

In order to produce plots/contours of Cp in Paraview:

	- import the injector solution, and select the latest tindx

	- SLICE the solution using Z-normal orientation, and then 
		set the Z origin to a very small number (0.00001)

	- using the CALCULATOR, enter the equation for Cp:
		(p-7100)/(0.5*(7100/(287.1*70.3))*(4*sqrt(1.4*287.1*70.3))^2)
		Note: p is a cell property, and requires selection from the list

	- add the CELL TO POINT DATA filter - the data can now be rescaled as
		desired, to indicate Cp over the flat plate surface!

	- use the CONTOUR tool to add Cp contours - Viti uses 13 entires, over
		the range of -0.1 to 0.5

############################
#### RUNNING A NEW CASE ####
############################

Tools have been included to generate new files for the 3D injector case, requiring
the running of a 2D run-up, interpolating a new inflow profile, extrapolating the
profile into 3D, and setting the generated files as inflow for the 3D simulation.
The interpolator was created quickly, and is likely to break when attempting to inter-
polate data from differently clustered/sized profiles. Use it at your own risk!
An easy solution is to simply make both 2D and 3D cases the same resolution, and 
just use to extrapolator.

A high resolution 2D `run-up' simulation should be run, and a lower resolution 3D
case. Alternatively, both can be the same resolution - this, however, will extend
the simulation time significantly.

Upon completion of both simulations, the 2D `run-up' should be sliced using the 
post.sh script, and a vertical slice of the 3D simulation should also be taken, to
provide cell data for the interpolator. 

Copy the .dat files into the interpolatorfolder and edit the python file as required. 
Executing the script will generate an interpolated 2D inflow profile. 

These profiles are then copied to the extrapolate folder. Edit the python file to 
supply the correct cell/block distribution, and execute it. This will produce the 
required number of 3D inflow .dat files (the extrapolator was built to produce 16 
block solutions - it should be simple to edit the file to produce more, if required.

Finally, the infiles should be set for use in the 3D simulation by editing the 
.config file, and the simulation can be executed!


----------------------------------------------------------------------------

