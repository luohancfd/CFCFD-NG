=======================================================
Simulation of stagnation streamline on FireII vehicle
=======================================================

.. Author: Daniel F. Potter
.. Date: Feb-2012

Description
-----------
This example uses Poshax3 to compute the
flow immediately behind the shock at
a flow condition which correspond to a
flight time of 1634 s
on the Fire II vehicle reentry trajectory.
The details of the calculation can be found in 
Chapter 5 of the author's thesis.

This example is found in the code collection at::

  $CFCFD/examples/poshax3/FireII/1634s

Files
-----

*User-built files:*

FireII.cfg
  Input file for Poshax3.

air-11sp-2T.inp
  Input file for describing the gas model.

Park93-s05-PEIIC.lua
  Input file for the chemical reaction model.

air-TV-TE.lua
  Input file for the energy exchange model.

run.sh
  A script to execute the example and plot results.

*Gas model input file:*

air-11sp-2T-HO.lua
  Gas model input file for direct use by Poshax3. This
  file is created using the ``gasfile`` program.

*Output file:*

output.data
  Output file from Poshax3 calculation. This contains flow properties
  as a function of distance downstream of the shock. *This file is only
  available after running the calculation.*

*Data files:*

marco_1634s\_??.txt
  Number densities and temperatures as calculated in Marco Panesi's PhD thesis.

*Results files:*

profiles.gplot
  A Gnuplot input file for plotting of results.

number_density_profiles.*
  Plot of number density profiles of key species.

temperature_profiles.*
  Plot of temperature profiles.

Running the simulation
----------------------

The steps to running and plotting results are contained in the script: ``run.sh``.
This can be executed at the command line, from this working directory as::

> ./run.sh

Results
-------
This example tests the ability of Poshax3 to correctly solve the post-shock
relaxation problem with both thermal and chemical nonequilibrium.
The output from the calculation are the flowproperties against distance behind
the shock.
A comparison of the temperature and number density profiles with the results
obtained by Marco Panesi in his PhD thesis are shown below.
The differences between the calculations can be attributed to variations in 
the thermodynamic modelling of key species such as N2.

.. figure:: temperature_profiles.png
   :align: center
   :scale: 150%

   Nonequilibrium temperatures behind a normal shock at Fire II t=1634s conditions.

.. figure:: number_density_profiles.png
   :align: center
   :scale: 150%

   Key species number densities behind a normal shock at Fire II t=1634s conditions.


