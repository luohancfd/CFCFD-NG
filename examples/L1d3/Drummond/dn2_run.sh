#! /bin/sh
# dn2_run.sh
# Sample script to exercise L1d on the Drummond tube 
# with Nitrogen driver.
#
# It is assumed that the L1d programs are correctly compiled and
# installed in a directory that is in the search path.
# Also, it is assumed that GNU-Plot is available as "gnuplot".

# Starting with the user's script dn2.py,
# create the input parameter file, dn2.Lp, and then 
# generate a starting solution file (dn2.L0) and 
# tube description (area) file (dn2.La).

l_script.py -f dn2
l_prep.exe -f dn2

# The simulation can now start at the initial (t=0) solution
# and integrate the gas-dynamic (and piston-dynamic) equations
# in time.
# The output files from this stage include:
# dn2.Ls  : the solution data for t>0
# dn2.hx  : history data for specific locations
# dn2.hc  : history data for specific (Lagrangian) gas cells
#
# Console output is captured in the file dn2.log which
# may be viewed during the simulation with the unix command tail -f.
#
# The Unix time command may be used to measure the CPU time required
# for the simulation.  The output from the time command seems to go
# to stderr and so is not caught in in the file dn2.log.

echo 
echo Beginning simulation...
echo Console output is being caught in file dn2.log.
time l1d.exe -f dn2 > dn2.log
echo Finished simulation.
echo 

date
