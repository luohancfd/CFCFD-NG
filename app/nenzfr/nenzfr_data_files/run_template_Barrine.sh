#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N $jobName
#PBS -q workq
#PBS -l place=scatter -l select=1:ncpus=1:NodeType=fast:mpiprocs=1 -l walltime=240:00:00 -A uq-SCRAMSPACE
#PBS -V

module load intel-mpi
module load intel-cc-11
module load python
echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "-------------------------------------------"
echo "Begin Shared job..."
date

#cd $HOME/Work/M10_LP_nominal/turb_grd100-1800-180blks_yclus1-002_6ms/neq_lam_throat/
#$HOME/Work/M10_LP_nominal/turb_grd100-1800-180blks_yclus1-002_6ms/neq_lam_throat/nenzfr.py --gas=air5species --T1=300.0 --p1=160.0e3 --Vs=2300.0 --pe=43.0e6 --cfile=contour-t4-m10.data --area=1581.165 --chem=neq > LOGFILE

nenzfr.py --gas=$gasName --T1=$T1 --p1=$p1 --Vs=$Vs --pe=$pe --chem=$chemModel --area=$areaRatio --job=$jobName --cfile=$contourFileName --gfile=$gridFileName --exitfile=$exitSliceFileName --nni=$nni --nnj=$nnj --nbi=$nbi --bx=$bx --by=$by --max-time=$max_time --max-step=$max_step > LOGFILE

echo "End Shared job."
date
