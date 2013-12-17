#! /bin/bash 
# run-optimiser.sh

echo "Start time: "; date
cd $PBS_O_WORKDIR
./optimise_T4_Mach_4_nozzle.py
echo "Finish time: "; date
