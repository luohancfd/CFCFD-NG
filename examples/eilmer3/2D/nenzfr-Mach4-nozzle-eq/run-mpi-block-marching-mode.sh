#! /bin/bash 
# run-mpi-block-marching-mode.sh
#
# Running nenzfr using the new multi-process block-marching mode 
# that has recently been implemented. Note the two new options
# needed when running in this mode (nbj and blockMarching variables
# that can be seen in the cfg file).

echo "Start time: "; date
nenzfr.py --config_file=nenzfr-Mach4-nozzle-eq-block-marching.cfg  > LOGFILE
echo "Finish time: "; date
