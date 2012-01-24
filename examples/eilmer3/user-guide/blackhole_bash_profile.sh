# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs
module load openmpi/1.2-gnu-4.1
module load swig
module load python/2.5.2
export PATH=${PATH}:${HOME}/bin:${HOME}/e3bin
export LUA_PATH=${HOME}/e3bin/?.lua
export LUA_CPATH=${HOME}/e3bin/?.so
unset USERNAME
