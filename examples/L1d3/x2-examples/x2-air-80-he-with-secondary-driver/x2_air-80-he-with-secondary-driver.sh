#! /bin/sh
# x2_air-80-he-with-secondary-driver.sh
# This is adapted from Dave Gildfind's current L1d script file. (Thanks Dave.)
#
# The python file coupled to this run script runs an L1d simulation
# of the X2 expansion tube with an 80%He/20%Ar primary driver, a pure He secondary driver
# and a nozzle. Air is the test case and the condition was just made up as an example.
#
# Chris James (c.james4@uq.edu.au) - 10/05/13
#
# It is assumed that the L1d programs are correctly compiled and
# installed in a directory that is in the search path.
#----------------------------------------------------------------

# Starting with the Python job file, x2_air-80-he-with-secondary-driver.py, 
# generate a starting solution file (x2_air-80-he-with-secondary-driver.L0) and 
# tube description (area) file (x2_air-80-he-with-secondary-driver.La).

l_script.py -f x2_air-80-he-with-secondary-driver
l1d.exe -f x2_air-80-he-with-secondary-driver -prep

# The simulation can now start at the initial (t=0) solution
# and integrate the gas-dynamic (and piston-dynamic) equations
# in time.
# The output files from this stage include:
# x2_air-80-he-with-secondary-driver.Ls  : the solution data for t>0
# x2_air-80-he-with-secondary-driver.hx  : history data for specific locations
# x2_air-80-he-with-secondary-driver.hc  : history data for specific (Lagrangian) gas cells
#
# Console output is captured in the file x2_air-80-he-with-secondary-driver.log which
# may be viewed during the simulation with the unix command tail -f.
#
# The Unix time command may be used to measure the CPU time required
# for the simulation. The output from the time command seems to go
# to stderr and so is not caught in in the file x2_air-80-he-with-secondary-driverlog.

echo 
echo Beginning simulation...
echo Console output is being caught in file x2_air-80-he-with-secondary-driver.log.
time l1d.exe -f x2_air-80-he-with-secondary-driver > x2_air-80-he-with-secondary-driver.log
echo Finished simulation.
echo 

date


