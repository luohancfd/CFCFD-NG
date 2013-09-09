#! /bin/bash 
# run.sh
#
e3prep.py --job=bump > LOGFILE_PREP
e3loadbalance.py --job=bump --number-of-procs=4
echo "Start time: "; date
mpirun -np 4 e3mpi.exe --job=bump --mpimap=bump.mpimap --run  > LOGFILE_RUN
echo "Finish time: "; date
e3post.py --job=bump --tindx=1 --add-mach --vtk-xml
