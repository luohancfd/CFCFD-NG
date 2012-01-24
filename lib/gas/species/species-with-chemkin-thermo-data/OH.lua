-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

OH = {}
OH.M = {
   value = 1.7007340e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
OH.gamma = {
   value = 1.3856e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
OH.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.992015430e+00, -2.401317520e-03, 4.617938410e-06, -3.881133330e-09, 1.364114700e-12, 3.615080560e+03, -1.039254580e-01, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.092887670e+00, 5.484297160e-04, 1.265052280e-07, -8.794615560e-11, 1.174123760e-14, 3.858657000e+03, 4.476696100e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
