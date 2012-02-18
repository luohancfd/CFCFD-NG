Poshax3
=======

Poshax3 is the latest installment of the Post-SHock relAXation solver first developed by Rowan 
Gollan.
The code solves for the spatial variation in flow properties behind a strong shock assuming 
inviscid and one-dimensional flow.
Source terms due to chemical reactions, thermal energy exchange and radiative cooling are able
to be considered.
A key feature of Poshax3 is the ability to treat these source terms in a loosely or fully
coupled manner.

Typical build and run procedure
-------------------------------
Poshax3 is built from source into an installation directory $HOME/e3bin/.  
A typical build procedure (using the default TARGET=for_gnu) might be::

  $ cd $HOME/cfcfd3/app/poshax3/build
  $ make install

Poshax3 calculations tend not to be very processor intensive, and therefore a parallelised 
version has not been developed. 
To see an example simulation, do:

  $ cd $HOME/cfcfd3/examples/poshax3/air-Mach-12.3
  $ ./run.sh
  
This shell script invokes the poshax3 calculation via:

  $ poshax3.x air-Mach-12.3.cfg
  
where air-Mach-12.3.cfg is the configuration script controlling the calculation.

Configuration script
--------------------
The parameters able to be set in a configuration script are listed in 
app/poshax3/source/poshax_spec.cfg.
Note the parameters must be located in their respective [sections].
A brief description of each is given here.

[controls]
dx                   : Initial or constant spatial step size
adaptive_dx          : Flag to indicate adaptive spatial stepping
final_x              : x location to terminate the calculation
plot_dx              : Spatial frequency of solution output 
dx_scale             : Scale the spatial step by this amount each step
output_file          : Name of the output file

[models]
gas_model_file       : Name of the gas model lua file
reaction_file        : Name of the reaction model lua file
energy_exchange_file : Name of the energy exchange model lua file
radiation_file       : Name of the radiation model lua file

[radiation]
rad_dx               : Spatial frequency of radiation intensity calculation
TS_dx                : Slab width for tangent-slab analysis
tube_width           : Width of gas to integrate over when computing intensities
fwhm_Ang             : Full-width half-maximum of the apparatus function in Angstroms 
x_EQ_spectra         : Location to calculate spectral coefficients and write to file
write_rad_level_pops : Flag to request writing of radiator electronic level populations to file
write_rad_emissions  : Flag to request writing of individual radiator total emission levels to file 
lambda_min           : Lower bound for spectral integration (can be a list)
lambda_max           : Upper bound for spectral integration (can be a list)
dx_smear             : Spatial smearing distance (i.e. shock speed x spectrometer exposure time)

[initial-conditions]
rho_inf              : Freestream density
p_inf                : Freestream pressure (if density not given)
T_inf                : Freestream temperature
u_inf                : Freestream velocity
M_inf                : Freestream Mach number (if velocity not given)
mass_f               : Freestream mass fractions
mole_f               : Freestream mole fractions (if mass fractions not given)

Examples
--------

.. toctree::
   :maxdepth: 1

   poshax3/air-Mach-12.3
   poshax3/FireII/1634s

