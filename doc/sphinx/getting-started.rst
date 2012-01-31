Getting the codes and preparing to run them
===========================================

The code repository
-------------------
The codes are available for download from a Mercurial repository.
To make a clone of the repository:

| $ cd $HOME
| $ hg clone https://cfcfdlocal@triton.pselab.uq.edu.au/cfcfd3-hg/cfcfd3-hg/ cfcfd3

This takes about 40 seconds on campus at UQ.  
It may take much longer, depending on your internet connection.

To see what's changed:

| $ cd cfcfd3
| $ hg incoming https://cfcfdlocal@triton.pselab.uq.edu.au/cfcfd3-hg/cfcfd3-hg/
| ...
| $ hg pull -u https://cfcfdlocal@triton.pselab.uq.edu.au/cfcfd3-hg/cfcfd3-hg/

Notes

#. You will need a password for any access.  Please ask.
#. You can read but not write with the "cfcfdlocal" username.
#. Some usernames (by negotiation) may push changesets back to the repository.


Your computational environment
------------------------------
The code collection comes as source code only so,
to use any of them, you will need to compile and install them.

To build and run the newer codes, you will need the following:

* a Unix-like system with GNU-make, C and C++ compilers
* popt (command-line parser) library and development files
* readline library (including the header files, libreadline5-dev on Ubuntu)
* Python + (with the numpy extension)
* SWIG
* Tcl/Tk + the BWidget library (to run the GUI program e3console.tcl)

We have been able to get the programs to build on Linux, MacOS-X 
(with a recent Xcode development environment) and Cygwin 1.7 (on MS-Windows).

On MS-Windows, install the full kit of Cygwin (Python, X-Windows and all)
and be careful not to have another Python installed outside of Cygwin.
The multiple installations of Python seem not to play well together.

Some other things that are useful:

* awk
* MetaPost (mpost) or, more recently, InkScape (for looking at and editing svg files)
* GNUplot
* Paraview or MayaVi or VisIt

To a basic Fedora 16 installation, you should add the following packages:

#. bzr
#. gcc
#. gcc-c++
#. m4
#. gcc-gfortran
#. swig
#. python-devel
#. readline-devel (for Lua)
#. popt-devel

To a basic Ubuntu 10.04 installation, you should add the following packages and their dependencies:

#. bzr
#. bzrtools
#. g++
#. m4
#. mpi-default-dev
#. gfortran
#. swig
#. python-dev
#. python-numpy
#. libreadline-dev
#. libpopt-dev
#. tk
#. bwidget
#. gnuplot


SSH access to the repository for developers
-------------------------------------------
Alternative access to the Mercurial repository for developers is possible via ssh.
You will need the password or your public key installed for any access.  Please ask.

| $ cd ~
| $ hg clone ssh://cfcfd3@triton/cfcfd3-hg cfcfd3
| $ cd cfcfd3/extern/
| $ hg clone ssh://cfcfd3@triton/cea2-hg cea2




