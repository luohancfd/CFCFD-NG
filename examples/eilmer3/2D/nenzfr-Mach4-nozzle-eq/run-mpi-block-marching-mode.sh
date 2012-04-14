#! /bin/bash 
# run-mpi-block-marching-mode.sh
#
# The new multi-process block-marching mode that has recently 
# been written into nenzfr.py. When running in this mode, ensure 
# that you have the nozzle.timing file that specifies the max_time 
# for each Eilmer3 run. 

echo "Start time: "; date
nenzfr.py --gas=air --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.376e6 --area=27.0 --cfile=contour-t4-m4.data \
          --nni=600 --nnj=40 --by=1.02 --nbi=60 --nbj=2 --block-marching > LOGFILE
echo "Finish time: "; date
