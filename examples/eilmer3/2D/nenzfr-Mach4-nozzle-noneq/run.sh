#! /bin/bash 
# run.sh

echo "Start time: "; date
nenzfr.py --config_file=nenzfr-Mach4-nozzle-neq.cfg > LOGFILE
echo "Finish time: "; date
