# PJ's bashrc for barrine

module load python
module load mercurial
module load intel-fc-13
module load intel-cc-13
module load intel-mpi
module load compiler/gcc-4.8.0
module load swig

export E3BIN=${HOME}/e3bin
export PATH=${PATH}:${E3BIN}
export LUA_PATH=${E3BIN}/?.lua
export LUA_CPATH=${E3BIN}/?.so
export PYTHONPATH=${PYTHONPATH}:${E3BIN}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${E3BIN}


