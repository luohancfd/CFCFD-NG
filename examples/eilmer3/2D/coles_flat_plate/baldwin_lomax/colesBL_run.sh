# colesBL_run.sh
date
e3prep.py --job=colesBL --do-svg > LOGFILE.prep
date
mpirun -np 4 e3mpi.exe --job=colesBL --run > LOGFILE.run
date
