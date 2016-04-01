Getting the codes and preparing to run them
===========================================

The code repository
-------------------
The codes are available for download from a public Mercurial repository.
To make a clone of the repository::

  $ cd $HOME
  $ hg clone https://source.eait.uq.edu.au/hg/cfcfd3 cfcfd3

You may provide the username "cfcfd-user@svn.itee.uq.edu.au" for authentication but
no password is required for reading.
Cloning the repository takes about 40 seconds on campus at UQ.  
It may take much longer, depending on your internet connection.

To see what's changed::

  $ cd cfcfd3
  $ hg incoming https://source.eait.uq.edu.au/hg/cfcfd3
  ...
  $ hg pull -u https://source.eait.uq.edu.au/hg/cfcfd3

Notes

#. You do not need a password for for read access.
#. You can read but not write with the "cfcfd-user" username.
#. Some usernames (cfcfd-dev@svn.itee.uq.edu.au) may push changesets back 
   to the repository.  You will need to negotiate a developer role for this access.
#. Some gas models depend on the NASA CEA code or the NIST REFPROP library.
   If you want to use these models (and there is no look-up-table equivalent
   already available) you will need to obtain these codes and place them 
   into the extern/ directory.  
   They are not included as part of our cfcfd3 repository but the cfcfd3 makefiles
   will be aware of them if they are sitting in the extern/ directory.


Licence
-------
CFCFD program collection is a set of flow simulation tools for compressible fluids.
Copyright (C) 1991-2016 Peter Jacobs, Rowan Gollan, Daniel Potter, 
Ingo Jahn, Anand Veeraragavan, Vince Wheatley, Daryl Bond, Chris James
and other members of the CFCFD group.

This collection is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program collection is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU-General-Public-License_
along with this program.  If not, see <http://www.gnu.org/licenses/>.

.. _GNU-General-Public-License: ./_static/gpl.txt


Your computational environment
------------------------------
The code collection comes as source code only so,
to use any of them, you will need to compile and install them.

To build and run the newer codes, you will need the following:

* a Unix-like system with GNU-make, C and C++ compilers
* popt (command-line parser) library and development files
* readline library (including the header files, libreadline-dev on Ubuntu)
* Python + the numpy, matplotlib and scipy extensions
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

#. mercurial
#. gcc
#. gcc-c++
#. m4
#. openmpi
#. openmpi-devel 
   (to use openmpi on Fedora, 
   :ref:`the module must be loaded <label-openmpi-fedora>`)
#. gcc-gfortran
#. libgfortran.i686, glibc-devel.i686 and libgcc.i686 
   (to compile the 32-bit CEA code on 64-bit Fedora)
#. swig
#. python-devel
#. numpy
#. python-matplotlib
#. scipy
#. readline-devel (for Lua)
#. popt-devel
#. sympy (to run the Method-of-Manufactured-Solutions test case for Eilmer3)

To a basic Ubuntu 10.04 (or any recent Debian derivative) installation, 
you should add the following packages and their dependencies:

#. mercurial
#. g++
#. m4
#. mpi-default-dev
#. mpi-default-bin
#. gfortran
#. gfortran-multilib (for compiling 32-bit CEA2 on a 64-bit system)
#. swig
#. python-dev
#. python-numpy
#. python-matplotlib
#. python-scipy
#. libreadline-dev
#. libpopt-dev
#. libncurses5-dev
#. tk
#. bwidget
#. gnuplot
#. tcl-dev (if you want to build IMOC)
#. python-sympy (to run the Method-of-Manufactured-Solutions test case for Eilmer3)

Compiler versions
-----------------
Since March 2013, we have started using some of the C++11 features 
such as range-based for loops and initializer expressions.
Because of this you will need a suitable C++ compiler.
For the GNU compiler collection, versions 4.6.3 and 4.8.0 are suitable.
Clang/LLVM versions 3.2 and later are also good.

Using the codes on MS-Windows
-----------------------------
The codes assemble most conveniently on a Linux/Unix-like environment.
They should also build and run within Cygwin (http://cygwin.com/), however,
it may be convenient to run a full linux installation within 
VirtualBox (https://www.virtualbox.org/), on your MS-Windows computer.

Using the codes on Apple OSX
----------------------------
The codes can be compiled and run on OSX as this is a Unix based OS.
The Xcode development environment (https://developer.apple.com/xcode/) 
should be downloaded and installed to provide Apple's versions of the 
GNU Compiler Collection, Python and the make utility, amongst other
development tools.
popt, readline, SWIG and Tcl/Tk can either be installed from source
or via a package manager such as MacPorts (http://www.macports.org/) or 
Fink (http://www.finkproject.org/).

Notes:

#. If possible, it is recommended to install these dependencies from source.
#. The required Python packages (numpy, scipy and matplotlib) are all available
   as pre-packaged binaries for OSX on sourceforge.net, although they can also
   be installed from source if necessary.
#. Ingo has had a good experience installing binary packages from MacPorts,
   the only subtly being the need to install swig and swig-python.

SSH access to the repository for developers
-------------------------------------------
Alternative access to the Mercurial repository for developers is possible via https.
You will need the password for the cfcfd-dev@svn.itee.uq.edu.au login.  Please ask.

::

  $ cd ~
  $ hg clone https://source.eait.uq.edu.au/hg/cfcfd3 cfcfd3
  $ cd cfcfd3/extern/
  $ hg clone https://source.eait.uq.edu.au/hg/cea2 cea2
  $ hg clone https://source.eait.uq.edu.au/hg/refprop refprop


