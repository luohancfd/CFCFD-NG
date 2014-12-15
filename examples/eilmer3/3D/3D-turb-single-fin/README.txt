# README.txt for 3D-turb-single-fin (tc2update.py)
# Test Case 2 - Three dimensional single fin utilising the k-omega turbulence model
# Author: Sam Stennett, Dec 2014

This is one of the test cases for Samuel Stennett's undergraduate thesis (2014).

The directory contains the images indicating the basic flow schematic, and the 
block arrangement used in the 3D simulation. 

#################################################################################
NOTE:
	As of the 15th Dec 2014, the case is still checker-boarding and results 
	in numerical instability. The case does not reach steady state, as flow
	conditions diverge to nan. Future work on this case is recommended, in 
	order to determine the origin of this phenomenon and prevent it.

#################################################################################

There are sufficient files in the directory to produce a 2D `run-up', interpolate
and extrapolate a 3D inflow profile, and simulate the complete 3D flow domain. 
As stated above, the current build of the case results in simulation instability
and checker-boarding. 

################################
#### RUNNING THE CASE AS IS ####
################################

To run the test case, enter the terminal and execute the following commands:

	e3prep.py --job=tc2update --verbosity=1

The verbosity flag allows for feedback during the prepping process, as the large
number of cells and blocks requires some time to prepare.
Next, the .config file must be edited to input the appropriate inflow profiles. 
If the .py script has NOT been altered, this can be done by executing the following:

	cp tc2update.config.backup tc2update.config

This will overwrite the current .config file with the version containing the new
boundary conditions. This can also be executed manually by editing the .config file
and replacing the REPLACE.dat entries with the appropriate .dat files. Currently, 
the geometry has 16 inflow boundaries (BLOCKS 0-15, WEST boundary), and can be set 
in ascending order repetitively:

	for B0: 	3Dinfilelower1.dat
	for B1: 	3Dinfilelower2.dat
	for B2: 	3Dinfileupper1.dat
	for B3: 	3Dinfileupper2.dat
	for B4: 	3Dinfilelower1.dat
	...		...
	for B15: 	3Dinfileupper2.dat

This is done as the blocks are defined vertically (0 to 3), and then step sideways
and repeat.
Finally, the case can be run by executing:

	./run_barrine.sh

The run_barrine.sh script is written for execution on the HPCU cluster `Barrine', 
or any system utilising PBS queue software and the MPI version of Eilmer. The expected 
runtime of the case is currently unknown, utilising 64 CPUs. 
If no MPI version is installed, or an insufficient number of CPUs, the case can 
instead be run using the command:

	e3shared.exe --job=tc2update --run

The runtime of the case executed in this manner will be very long (months?) on a 
typical PC or laptop. It is highly recommended that `Barrine' access is granted in
order to run the 3D simulations, due to the large number of cells and use of the 
k-omega turbulence model.

############################
#### RUNNING A NEW CASE ####
############################

Tools have been included to generate new files for the 3D single fin case, requiring
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

