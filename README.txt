The Compressible-Flow CFD code collection
=========================================
:Author: Peter_Jacobs and Rowan_Gollan
:Email: peterj@mech.uq.edu.au 


The codes
---------
.The main codes are found in the app subdirectory.
- *bos*: Background-oriented schlieren (current but needs maintenance; new numpy)
- *elmer*: The first version of our time-dependent 3D flow simulation code.
  C version built for LES simulations. (obsolete)
- *elmer2*: Second version using libgas2 thermochemsitry and written in C++.
  (obsolete)
- *eilmer3*: The current compresible-flow Navier-Stokes simulation code
  for 2D and 3D flows. 
  This has the best bits of mbcns2 and elmer2, together with a rebuilt
  thermochemistry module.
  We have changed the name to Eilmer, the correct spelling of the name of the
  patron monk of CFD and aeronautics.
  See the wikipedia article http://en.wikipedia.org/wiki/Eilmer_of_Malmesbury
- *imoc*: interactive method-of-characteristics code, written in C+Tcl.
- *L1d*: Lagrangian simulation of quasi-one-dimensional flow. (obsolete)
- *L1d2*: second version of this code build with libgas2 thermochemistry.
  Written in C++. (current)
- *mb_cns*: Multiple-block compressible Navier-Stokes flow simulator.
  This is the code corresponding to the 1996 report. (obsolete)
- *mbcns2*: second version using libgas2 and written in C++. (obsolete)
- *newt*: Newtonian flow approximation for hypersonics flows.
- *octvce*: Joseph's fast Euler code with hierarchical adaptive grid for 3D flows.
- *poshax*: one-dimensional flow with complex thermochemistry.
- *sdas*: firmware for simple data acquisition system for student Zuni projects.
- *sf3d_v4.1*: Ian Johnston's shock-fitting code for 3D flow.
- *sm_3d_plus*: space-marching code for 3D flow.  
  (obsolete now that block-sequencing is part of eilmer3)
- *tds*: tunnel data server, software for both the server and the client


Getting the simulation codes
----------------------------
The top-level repository doesn't have many files in itself.
Most of the interesting stuff is in the sub-repositories:
app, doc, examples and lib.
One may get a complete working copy via http-based checkout
(with read-only access).

[source,shell]
cd $HOME
svn checkout --username cfcfd http://triton.pselab.uq.edu.au/repos/cfcfd2 cfcfd2
cd cfcfd2
svn checkout http://triton.pselab.uq.edu.au/repos-cfcfd2/app app
svn checkout http://triton.pselab.uq.edu.au/repos-cfcfd2/lib lib
svn checkout http://triton.pselab.uq.edu.au/repos-cfcfd2/examples examples
svn checkout http://triton.pselab.uq.edu.au/repos-cfcfd2/doc doc

.Notes
- When prompted, the password is: cfd-socks
- There is a hyphen between repos and cfcfd2 for the second through fifth
checkouts.


Getting the simulation codes -- Developers
------------------------------------------
Developers should get the most up-to-date code using bazaar-ng.
To get (read-only) access via http, use the command:

[source,shell]
cd $HOME
bzr branch http://cfcfdlocal@triton.pselab.uq.edu.au/bzr-cfcfd2 cfcfd2

.Notes
- You will need a password for bzr access.  Please ask.
- To push code changes, you will need access via ssh.


Your build environment
----------------------
The code collection comes as source code only so,
to use any of them, you will need to compile and install them.

.To build and run the newer codes, you will need the following:
. a Unix-like system with GNU-make, C and C++ compilers
. popt (command-line parser) library and development files
. readline library (including the header files, libreadline5-dev on Ubuntu)
. Python + (with the numpy extension)
. SWIG
. Tcl/Tk + the BWidget library (to run the GUI program e3console.tcl)

We have been able to get the programs to build on Linux, MacOS-X 
(with a recent Xcode development environment) and Cygwin 1.7 (on MS-Windows).

On MS-Windows, install the full kit of Cygwin (Python, X-Windows and all)
and be careful not to have another Python installed outside of Cygwin.
The multiple installations of Python seem not to play well together.

.Some other things that are useful:
. awk
. MetaPost (mpost) or, more recently, InkScape (for looking at and editing svg files)
. GNUplot
. Paraview or MayaVi or VisIt

.To a basic Fedora 16 installation, you should add the following packages:
. bzr
. gcc
. gcc-c++
. m4
. gcc-gfortran
. swig
. python-devel
. readline-devel (for Lua)
. popt-devel

.To a basic Ubuntu 10.04 installation, you should add the following packages and their dependencies:
. bzr
. bzrtools
. g++
. m4
. mpi-default-dev
. gfortran
. swig
. python-dev
. python-numpy
. libreadline-dev
. libpopt-dev
. tk
. bwidget
. gnuplot


Typical build and run procedures
--------------------------------

mbcns2 (2D flow simulation)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
[source,shell]
cd $HOME/cfcfd2/app/mbcns2/build
make
make install
make clean
cd $HOME/cfcfd2/lib/plot/build
make
make install
make clean

You may need to add the installation directory to your system's 
search path and to Lua's search path.
On a recent Linux system, this could be done by adding the lines

[source,shell]
export PATH=${PATH}:${HOME}/cfd_bin
export LUA_PATH=${HOME}/cfd_bin/?.lua
export LUA_CPATH=${HOME}/cfd_bin/?.so

to the .bash_profile file in your home directory.

Then, try out the cone20 example.

[source,shell]
mkdir $HOME/work; cd $HOME/work
mkdir cone20; cd cone20
cp $HOME/cfcfd2/examples/mbcns2/cone20/* .
./cone20_run.sh
./cone20_plot.sh

This should generate a few postscript figures of the starting flow
about a sharp 20-degree cone.

Elmer2 
~~~~~~
This 3D-flow code is built similarly.  
For running on the blackhole cluster:

[source,shell]
cd $HOME/cfcfd2/app/elmer2/build
module load openmpi/1.2-gnu-4.1
module load swig
module load lua
make TARGET=for_openmpi WITH_MPI=1 install

To run the flow over a cylinder simulation with dissociating nitrogen:

[source,shell]
mkdir $HOME/work/elmer2; cd $HOME/work/elmer2
mkdir tim1; cd tim1
cp $HOME/cfcfd2/examples/elmer2/tim1/* .
./prepare_thermo_files.sh
qsub cyl_run.sh

Eilmer3
~~~~~~~
The new 2D/3D code Eilmer3 is built much more like mbcns2 but
has a new installation directory $HOME/e3bin/.  
A typical build procedure (using the default TARGET=for_gnu) might be:

[source,shell]
cd $HOME/cfcfd2/app/eilmer3/build
make install
make clean

Or, if you want the MPI version of the code built as well:

[source,shell]
cd $HOME/cfcfd2/app/eilmer3/build
make TARGET=for_openmpi install
make clean


You may need to add the installation directory to your system's 
search path and to Lua's search path.
On a recent Linux system, this could be done by adding the lines

[source,sh]
export PATH=${PATH}:${HOME}/e3bin
export LUA_PATH=${HOME}/e3bin/?.lua
export LUA_CPATH=${HOME}/e3bin/?.so

to the .bash_profile or .profile file in your home directory.
(Note that it is not necessary to have the LUA_CPATH variable set 
unless you want to access the Lua gas module from within 
the user-defined (Lua) functions.
You don't need it to run Eilmer3 otherwise.)

For running on the blackhole cluster, also add the following:

[source,sh]
module load openmpi/1.2-gnu-4.1
module load swig

IMPORTANT: Do NOT load lua as was done in the example for Elmer2.

Then, try out the cone20-simple example.

[source,sh]
mkdir $HOME/work; cd $HOME/work; mkdir 2D; cd 2D
mkdir cone20-simple; cd cone20-simple
cp $HOME/cfcfd2/examples/eilmer3/2D/cone20-simple/* .
./cone20_run.sh  # exercise the shared-memory version of the code
       or
./cone20_run_mpi.sh  # exercise the MPI version of the code

This should generate a postscript figure of the drag coefficient history
about a sharp 20-degree cone and also put the VTK data file into the plot/
subdirectory.
It is not really necessary to make all of the subdirectories as shown above,
however, that arrangement reflects the directory tree that PJ uses.
If you want him to come and look at your simulation files when things go wrong,
use the same.
If not, use whatever hierarchy you like.


Other Notes
-----------
.On Xserver for Linux (especially Ubuntu)
If Paraview crashes on exporting a bitmap image, try adding the line

[source,shell]
Option "AIGLX" "false"

to the Section "ServerLayout" in /etc/X11/xorg.conf

.To use Paraview 3.6.1 on Ubuntu 9.04 or later,
it seems that we need to customize the look of the desktop 
by turning off the Visual Effects. 
This setting can be found in the System->Preferences->Appearance menu.

.To get Paraview Screenshot to behave,
uncheck "Use Offscreen Rendering for Screenshots" button
in the Edit->Settings ("Options") dialog.
You will find the checkbutton under "Render View"->General.
