#! /bin/bash 
# run.sh

echo "Start time: "; date
nenzfr.py --gas=air5species --chem=neq --T1=300 --p1=125.0e3 --Vs=2414 --pe=34.376e6 --area=27.0 --cfile=contour-t4-m4.data \
          --nni=600 --nnj=40 --by=1.02 --nbi=60 --max-time=1.0e-3 > LOGFILE
echo "Finish time: "; date
