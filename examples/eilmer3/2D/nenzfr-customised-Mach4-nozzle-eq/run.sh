#! /bin/bash 
# run.sh

echo "Start time: "; date
nenzfr.py --config_file=nenzfr-customised-Mach4-nozzle-eq.cfg > LOGFILE
echo "Finish time: "; date
