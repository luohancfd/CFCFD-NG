#! /bin/bash 
# run-nenzfr-in-mpi-block-marching-mode.sh

echo "Start time: "; date
nenzfr.py --gas="air" --T1=300 --p1=131.0e3 --Vs=1644.0 --pe=8.32e6 --chem="eq" \
          --area=129.72 --cfile=Bezier-control-pts-t4-m7.data --nni=300 --nnj=40 --by=1.002 --bx=1.5 \
          --nbi=30 --nbj=4 --block-marching --max-time=2.0e-3 --max-step=1000000 --BLTrans=0.00
echo "Finish time: "; date

