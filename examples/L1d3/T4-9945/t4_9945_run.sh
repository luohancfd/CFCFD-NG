#! /bin/sh
# t4_9945_run.sh

# Starting with the Python job file, t4_9945.py, 
# generate a starting solution file (t4_9945.L0) and 
# tube description (area) file (t4_9945.La).

l_script.py -f t4_9945
l1d.exe -f t4_9945 -prep

# The simulation can now start at the initial (t=0) solution
# and integrate the gas-dynamic (and piston-dynamic) equations
# in time.
# The output files from this stage include:
# t4_9945.Ls  : the solution data for t>0
# t4_9945.hx  : history data for specific locations
# t4_9945.hc  : history data for specific (Lagrangian) gas cells
#
# Console output is captured in the file t4_9945.log which
# may be viewed during the simulation with the unix command tail -f.
#
# The Unix time command may be used to measure the CPU time required
# for the simulation.  The output from the time command seems to go
# to stderr and so is not caught in in the file t4_9945.log.

echo 
echo Beginning simulation...
echo Console output is being caught in file t4_9945.log.
time l1d.exe -f t4_9945 > t4_9945.log
echo Finished simulation.
echo 

date
