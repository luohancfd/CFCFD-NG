#! /bin/bash 
# run-mpi-block-marching-mode.sh
#
# Running nenzfr using the new multi-process block-marching mode 
# that has recently been implemented. Note the two new options
# needed when running in this mode --nbj and --block-marching.

echo "Start time: "; date
nenzfr.py --gas=air --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.376e6 --area=27.0 --cfile=Bezier-control-pts-t4-m4.data \
          --nni=600 --nnj=40 --by=1.02 --nbi=60 --nbj=2 --block-marching --max-time=1.0e-3 > LOGFILE
echo "Finish time: "; date
