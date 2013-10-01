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
To see an example simulation, do::

  $ cd $HOME/cfcfd3/examples/poshax3/air-Mach-12.3
  $ ./run.sh
  
This shell script invokes the poshax3 calculation via::

  $ poshax3.x air-Mach-12.3.cfg
  
where air-Mach-12.3.cfg is the configuration script controlling the calculation.

Configuration script
--------------------
The parameters able to be set in a configuration script are listed in 
app/poshax3/source/poshax_spec.cfg.
Note the parameters must be located in their respective [sections].
A brief description of each is given here.

::

  [controls]
  source_term_coupling = string(default=loose)
  dx                   = float(default=1.0e-9)
  adaptive_dx          = boolean(default=true)
  final_x              = float(default=1.0e-3)
  plot_dx              = float(default=final_x*1.0e-3)
  dx_scale             = float(default=1.0)
  output_file          = float(default=poshax_output.data)
  apply_udpedx         = boolean(default=false)

  [models]
  gas_model_file       = string(default=None)
  reaction_file        = string(default=None)
  energy_exchange_file = string(default=None)
  radiation_file       = string(default=None)

  [radiation]
  rad_dx               = float(default=0.0)
  TS_dx                = float(default=0.0)
  path_length          = float(default=0.0)
  fwhm_Ang             = float(default=0.0)
  write_rad_level_pops = boolean(default=false)
  write_rad_emissions  = boolean(default=false)
  write_spectra        = boolean(default=false)
  lambda_min           = float_list(default=None)
  lambda_max           = float_list(default=None)
  dx_smear             = float_list(default=None)

  [initial-conditions]
  rho_inf      = float(default=None)
  p_inf        = float(default=None)
  T_inf        = float(default=None)
  u_inf        = float(default=None)
  M_inf        = float(default=None)
  mass_f       = float_list(default=None)
  mole_f       = float_list(default=None)

Examples
--------

.. toctree::
   :maxdepth: 1

   poshax3/air-Mach-12.3
   poshax3/thivet_et_al
   poshax3/FireII-1634s
   poshax3/Glass-Liu

