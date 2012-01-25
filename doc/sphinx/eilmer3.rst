Eilmer3
=======

Eilmer3 is our principal simulation code for 2D and 3D gas dynamics.
It is a research and education code, suitable for the exploration of
flows where the bounding geometry is not too complex.

Typical build and run procedure
-------------------------------
The new 2D/3D code Eilmer3 is built from source into an installation directory $HOME/e3bin/.  
A typical build procedure (using the default TARGET=for_gnu) might be:

| cd $HOME/cfcfd2/app/eilmer3/build
| make install
| make clean

Or, if you want the MPI version of the code built as well:

| cd $HOME/cfcfd2/app/eilmer3/build
| make TARGET=for_openmpi install
| make clean

You may need to add the installation directory to your system's 
search path and to Lua's search path.
On a recent Linux system, this could be done by adding the lines

| export PATH=${PATH}:${HOME}/e3bin
| export LUA_PATH=${HOME}/e3bin/?.lua
| export LUA_CPATH=${HOME}/e3bin/?.so

to the .bash_profile or .bashrc file in your home directory.
Note that it is not necessary to have the LUA_CPATH variable set 
unless you want to access the Lua gas module from within 
the user-defined (Lua) functions.
You don't need it to run Eilmer3 otherwise.

For running on some managed computers, such as our blackhole cluster, also add the following:

| module load openmpi/1.2-gnu-4.1
| module load swig

IMPORTANT: Do NOT load lua as was done in the example for Elmer2.

Then, try out the cone20-simple example.

| mkdir $HOME/work; cd $HOME/work; mkdir 2D; cd 2D
| mkdir cone20-simple; cd cone20-simple
| cp $HOME/cfcfd2/examples/eilmer3/2D/cone20-simple/* .
| ./cone20_run.sh  # exercise the shared-memory version of the code
|        or
| ./cone20_run_mpi.sh  # exercise the MPI version of the code

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
On Xserver for Linux (especially Ubuntu):

* If Paraview crashes on exporting a bitmap image, try adding the line
  
  Option "AIGLX" "false"

  to the Section "ServerLayout" in /etc/X11/xorg.conf

* To use Paraview 3.6.1 on Ubuntu 9.04 or later,
  it seems that we need to customize the look of the desktop 
  by turning off the Visual Effects. 
  This setting can be found in the System->Preferences->Appearance menu.

* To get Paraview Screenshot to behave,
  uncheck "Use Offscreen Rendering for Screenshots" button
  in the Edit->Settings ("Options") dialog.
  You will find the checkbutton under "Render View"->General.

User Guide (PDF)
----------------




