=======================================================
Simulation of stagnation streamline on FireII vehicle
=======================================================

.. Author: Daniel F. Potter
.. Date: Oct-2013

Description
-----------
This example uses Poshax3 to reproduce the ionisation 
fraction measurements made by Glass and Liu (JFM 1978
vol. 84 part1 pp. 55-77) in the UTIAS shock tube with 
an Argon test gas. The comparison is for the Ms=16.5 
case.

This example is found in the code collection at::

  $CFCFD/examples/poshax3/Argon/Glass-Liu/M16.5
  
where there are cases both with and without radiation coupling::

  $CFCFD/examples/poshax3/Argon/Glass-Liu/M16.5/with_radiation
  $CFCFD/examples/poshax3/Argon/Glass-Liu/M16.5/without_radiation

The files list below is the case with radiation coupling.

Files
-----

*User-built files:*

argon.cfg
  Input file for Poshax3.

argon-3sp-2T.inp
  Input file for describing the gas model.

Hoffert67-Ar-2T.lua
  Input file for the chemical reaction model.

Ar-TE.lua
  Input file for the energy exchange model.
  
argon-radiators-NIST-TB.py
  Input file for describing the radiation model.

run.sh
  A script to execute the example and plot results.

*Gas model input file:*

argon-3sp-2T.lua
  Gas model input file for direct use by Poshax3. This
  file is created using the ``gasfile`` program.
  
*Radiation model input file:*

rad-model.lua
  Radiation model input file for direct use by Poshax3. This
  file is created using the ``radmodel.py`` program.

*Output file:*

output.data
  Output file from Poshax3 calculation. This contains flow properties
  as a function of distance downstream of the shock. *This file is only
  available after running the calculation.*

*Data files:*

glass_liu_ne_M_16.5.txt
  Ionisation fraction from Figure 7 of Glass and Liu (1978).

*Results files:*

profiles.gplot
  A Gnuplot input file for plotting of results.

ionization_fraction_profile.*
  Plot of the ionisation fraction (electron mole fraction) profile.

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
relaxation problem with both thermochemical nonequilibrium and radiation coupling.
The output from the calculation are the flow properties against distance behind 
the shock.
The temperature profile for the calculation with radiation coupling is shown below.
The heavy particle translation temperature (black) is much larger than the electron 
temperature immediately behind the shock, as is assumed the translational excitation
of free electrons is much slower than that for heavy particles.
The two temperature equilibriate approximately 1.75cm behind the shock, and subsequently 
drop gradually with distance due to radiation cooling.

.. figure:: GL_temperature_profiles.png
   :align: center
   :scale: 150%

   Heavy particle translation and electron temperatures behind the shock.

Comparisons of the ionisation fraction profiles with the results obtained by Glass
and Liu in the UTIAS shock tube via two-wavelength interferometry are presented below.
The observed relaxation length of 1.8cm is only slightly larger than that predicted by
the calculations, however the peak ionisation fraction is slightly underestimated.
This is possibly due to deviations of the fill pressure and shock speed from the nominal
values in the experiment.
A clear difference between the calculations with and without radiation coupling is the 
slope of the ionisation fraction after the thermochemical equilibration point.
Radiation coupling is seen to result in a constant drop in the ionisation fraction, while
without radiation coupling the ionisation fraction maintains a steady equilibrium value.
The similarity of the observed and calculated slope in this region indicated the 
radiation model is accurately predicting the rate of cooling.  
In these calculations radiation at wavelengths less than 200nm was considered to be 
optically thick, and that at higher wavelengths optically thin.
This is a reasonable assumption as VUV radiation is strongly self absorbed.  

In the first figure, radiation coupling is considered and 

.. figure:: GL_alpha_with_radiation.png
   :align: center
   :scale: 150%

   Comparison of calculated (with radiation coupling) and measured ionisation fraction
   behind a normal shock in argon at Mach 16.5.

.. figure:: GL_alpha_without_radiation.png
   :align: center
   :scale: 150%

   Comparison of calculated (without radiation coupling) and measured ionisation fraction
   behind a normal shock in argon at Mach 16.5.
