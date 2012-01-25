Getting Started
===============

Getting the current code collection
-----------------------------------
The codes are available for download via a Mercurial repository.

| $ cd ~
| $ hg clone ssh://cfcfd3@triton/cfcfd3-hg cfcfd3
| $ cd cfcfd3/extern/
| $ hg clone ssh://cfcfd3@triton/cea2-hg cea2

Notes

#. You will need a password for any access.  Please ask.
#. To push code changes, you will need access via ssh.


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


Getting the cfcfd2 (old) codes
------------------------------
The older cfcfd2 code collection is still available using bazaar-ng.
To get (read-only) access via http, use the command:

| cd $HOME
| bzr branch http://cfcfdlocal@triton.pselab.uq.edu.au/bzr-cfcfd2 cfcfd2

Notes

#. You will need a password for bzr access.  Please ask.
#. To push code changes, you will need access via ssh.


