======================================
Cheat Sheet for Using Eilmer on Magnus
======================================

Logging in
----------

The hostname for the supercompute cluster is ``magnus.pawsey.org.au``.

Login is via ssh. From a linux terminal, for example::

  > ssh username@magnus.pawsey.org.au


Setting your computational environment
--------------------------------------

For INTEL compiler:

Add these lines to the end of your ``.bashrc``::

  module swap PrgEnv-cray PrgEnv-intel
  module swap craype-haswell craype-sandybridge
  module load swig
  module load mercurial
  
  export PATH=${PATH}:${HOME}/e3bin
  export LUA_PATH=${HOME}/e3bin/?.lua
  export LUA_CPATH=${HOME}/e3bin/?.so
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/e3bin

For GNU compiler:

Add these lines to the end of your ``.bashrc``::

  module swap PrgEnv-cray PrgEnv-gnu
  module swap craype-haswell craype-sandybridge
  module load swig
  module load mercurial
  module load scipy
  
  export PATH=${PATH}:${HOME}/e3bin
  export LUA_PATH=${HOME}/e3bin/?.lua
  export LUA_CPATH=${HOME}/e3bin/?.so
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/e3bin
  
GNU Compiler is compatible with module scipy and is selected
when you want to use module scipy, otherwise Intel compiler is recommended.    

This sets your environment ready for compiling and running.
It's easiest to log out and back in to have the changes take effect.

Compiling eilmer
----------------

Get a fresh copy of the repository::

  > hg clone https://source.eait.uq.edu.au/hg/cfcfd3 cfcfd3

Username: cfcfd-user@svn.itee.uq.edu.au
Password: hyper-socks

Copy the ``cea2`` source directory to ``cfcfd3/app/extern/``.

Compile eilmer::

  > cd cfcfd3/app/eilmer3/build
  > make TARGET=for_cray_intel install
  or
  > make TARGET=for_cray_gnu install
  
Submitting a job
----------------

You should run your jobs from the ``/group`` disk, so copy your working files over to that partition.
(It's ok to leave the source code in your home area).

Your ``/group`` directory is of the form::

  /group/"project-no"/"username"

For example, my ``/group`` directory is::
  
  /group/fy8/kqin
  
The actual job submission would look something like::

  > cd /group/fy8/kqin
  > cd sim_dir
  --- prepare submit script ---
  > qsub submit-script

An example submit script is::

  #!/bin/bash --login
  
  #SBATCH --nodes=2
  #SBATCH --time=24:00:00
  #SBATCH --account=fy8
  #SBATCH --output=hostname.out
  #SBATCH --error=hostname.err
  
  echo "Start MPI job...."
  date
  aprun -n 48 e3mpi.exe --job=jobname --run > LOGFILE
  echo "End MPI job."
  date

Note that on Magnus, a node has 24 processors. If you request more than one node, then
you need to request in multiples of 24, even if you don't use all of those processors.
The account will be charged though as if you used all processors requested.
This is not really a problem using eilmer: you can use the SuperBlock facility
to split one large block into many smaller blocks and the ``e3loadbalance`` utility
along with the ``--mpimap`` option to run the code on a desired number of
processors (which will typically be a multiple of 24).

Account usage
-------------
The available compute hours in each quarter, and how many of those
hours have been used are available via the ``pawseyAccountBalance`` command:

  > pawseyAccountBalance -users

More information
----------------
A User Guide for Magnus::

  https://portal.pawsey.org.au/docs/Supercomputers/Magnus/Magnus_User_Guide




  
