#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N $caseName
#PBS -q workq
#PBS -l place=scatter -l select=1:ncpus=8:NodeType=medium:mpiprocs=8 -l walltime=50:00:00 -A uq-SCRAMSPACE
#PBS -V

# With the exception of the above PBS options, this file should not need to be modified.
# Luke Doherty
# 16-04-2012

module load intel-mpi
module load intel-cc-11
module load python
echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "-------------------------------------------"
echo "Begin Shared job..."
date

cd $PBS_O_WORKDIR
$HOME/e3bin/nenzfr.py --gas=$gasName --T1=$T1 --p1=$p1 --Vs=$Vs --pe=$pe --chem=$chemModel --area=$areaRatio --job=$jobName --cfile=$contourFileName --gfile=$gridFileName --exitfile=$exitSliceFileName $blockMarching --nni=$nni --nnj=$nnj --nbi=$nbi --nbj=$nbj --bx=$bx --by=$by --max-time=$max_time --max-step=$max_step --Twall=$Tw --BLTrans=$BLTrans --TurbVisRatio=$TurbVisRatio --TurbIntensity=$TurbInten --CoreRadiusFraction=$CoreRadiusFraction > LOGFILE

echo "End Shared job."
date
# As we leave the job, make sure that we leave no processes behind.
# (The following incantation is from Gerald Hartig.)
for i in $(cat $PBS_NODEFILE | grep -v 'hostname' | sort -u); do
    ssh $i pkill -u 'whoami'
done
killall -u 'whoami' e3mpi.exe
