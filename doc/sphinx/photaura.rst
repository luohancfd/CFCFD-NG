Photaura
========

Photaura is a spectral radiation model that can be used either to perform 
spectral analyses via a standalone python library, or used to compute the 
spectral coefficients required for radiation transport calculations 
with the Eilmer3 and Poshax3 codes.
The name Photaura comes from the Greek word photos (light) and the Latin
word aura (wind); an appropriate name for a code that deals with plasmas.

Typical build and test procedure
--------------------------------
The parts of the Photaura code that are used by Eilmer3 and Poshax3 are
automatically compiled and linked by the these codes own makefiles.
To work with the standalone library, however, the user needs to do
'make install' in the radiation build area:

  $ cd $HOME/cfcfd3/lib/radiation/build
  $ make install

The installation can then be verified by running an automated test script:

  $ cd $HOME/cfcfd3/lib/radiation/build
  $ make test

Note that a number of the python tools provided in the radiation library 
make use of the cea2 interface located in lib/cfpylib/gasdyn/cea2_gas.py.
The successful use of this interface requires the cea2 source file be 
be placed in cfcfd3/extern/cea2 when building the code
(see http://www.grc.nasa.gov/WWW/CEAWeb/).
Also, to enable plotting to the screen the matplotlib python module
(which provides pylab) should also be installed 
(see http://matplotlib.sourceforge.net/).

Python input script
-------------------
See the file cfcfd3/lib/radiation/test/air-radiators.py for an example
input script.
More information coming soon.


Examples
--------

.. toctree::
   :maxdepth: 1

   photaura/VKI-minitorch/README

