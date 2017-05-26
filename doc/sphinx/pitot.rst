PITOT
=======

PITOT is the group's equilibrium gas full facility simulation tool for rapid characterisation and simulation of generally free piston driven expansion tube flow conditions. It makes use of isentropic and compressible flow relations to rapidly simulation the facility, while also including the ability to analytically account for some situations where expansion tube facilities are known to depart from ideal shock tube theory.

PITOT was written to be a virtual impulse facility, and simulations are therefore configured like an experimenter would configure a real experiment. It uses facility fill condition as inputs and then the code runs through the flow processes in a state-to-state manner, analogous to how the different sections of the facility would operate in the real experiment. The code was written this way to create a simple and intuitive tool for trying to understand a facility and how different parameters affect flow conditions. PITOT can also be easily scripted to perform parametric studies and sensitivity analyses, and tools are provided with the code to do this.

PITOT is written in the Python programming language and it is generally ran by parsing a configuration file to the overarching pitot.py program, which is discussed below. PITOT can also be scripted from inside Python by parsing a configuration dictionary directly to the run_pitot function which can be found in the overarching program pitot.py.

Publications about the code
-------------------------------

Several different papers have been written about the code, discussing what it does and how it works.  

The earliest paper was written in 2013:

James, C.M., Gildfind, D.G., Morgan, R.G., Jacobs ,PA and Zander, F. (2013) Designing and simulating high enthalpy expansion tube conditions. In: Proceedings of the 2013 Asia-Pacific International Symposium on Aerospace Technology, Takamatsu, Japan, (1-10). 20-22 November 2013.

Section 3 of this paper from 2015 also provides a more up to date discussion of the code:

James, C.M., Gildfind, D.G., Morgan, R.G., Lewis, S.W., Fahy, E.J., and McIntyre, T.J. (2015) On the current limits of simulating gas giant entry flows in an expansion tube. In: 20th AIAA International Space Planes and Hypersonic Systems and Technologies Conference. AIAA International Space Planes and Hypersonic Systems and Technologies Conference, Glasgow, Scotland, (1-26). 6 - 9 July 2015. doi:10.2514/6.2015-3501

However, by far the most comprehensive discussion of the code is presented in this paper from 2017:

James, C.M., Gildfind, D.G., Lewis, S.W., Morgan, R.G., and Zander, F. Implementation of a state-to-state analytical framework for the calculation of expansion tube flow properties. Submitted to Shock Waves, 2017.

Typical build and run procedure
-------------------------------
PITOT is built from source into a default installation directory at $HOME/e3bin/.  
A typical build procedure (using the default TARGET=for_gnu) might be::

  $ cd $HOME/cfcfd3/app/pitot
  $ make install

It should be noted that cea2 files must be in the cfcfd3/extern/cea2 directory when Pitot is compiled for the build to work correctly and allow PITOT to be ran in the normal equilibrium gas mode. The instructions for obtaining cea2 locally, as well as the dependencies required to run Pitot and other cfcfd programs, can be found `here <getting-started.html>`_. The uncompiled version of cea2 can also be obtained directly from NASA `here <https://www.grc.nasa.gov/WWW/CEAWeb/ceaguiDownload-unix.htm>`_ by downloading the file `CEA+Fortran.tar.Z'. These files must then be placed in the cfcfd3/extern/cea2 directory.

Equilibrium gas pitot simulations generally take a few minutes to run.
Several example simulations for various situations exist in the examples directory of the CFCFD repository (cfcfd3/examples/).
A basic example simulation can be ran with the commands below::

  $ cd $HOME/cfcfd3/examples/pitot/high-speed-air-theory
  $ ./high-speed-air-theory.sh
  
This shell script runs the pitot calculation by telling the program where the simulation cfg file is::

  $ pitot.py --config_file=high-speed-air-theory.cfg
  
Configuration script
--------------------
A comprehensive set of configuration options for PITOT is not yet compiled. 
Currently, information about configuring the code for various situations can be found in the various examples found in the cfcfd3/examples/ directory.

The configuration file is loaded into Python as a dictionary, so it must conform to Python syntax. Inside the configuration file, variables can be placed in any order, and Python's comment character (#) can be used to add text comments to the configuration file.
