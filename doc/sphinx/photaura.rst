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

In order to run the test programs, the test target needs to be requested:

  $ cd $HOME/cfcfd3/lib/radiation/build
  $ make test

Note that for best performance of the test programs the cea2 source code
should be placed in cfcfd3/extern/cea2.
The 'make test' command will then build the cea2 executable so that 
it can be used via the cea2_gas.py interface.
Also, to enable plotting to the screen the matplotlib python module
(which provides pylab) should also be installed 
(see http://matplotlib.sourceforge.net/).

To run the test program after running 'make test':

  $ cd $HOME/cfcfd3/lib/radiation/test
  $ ./run_test.sh

This will run some simple calculations with the radpy module
that tests the basic features of a Photaura spectral model.

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

