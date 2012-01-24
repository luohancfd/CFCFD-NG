-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

H2 = {}
H2.M = {
   value = 2.0158800e-03,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
H2.gamma = {
   value = 1.4049e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
H2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.344331120e+00, 7.980520750e-03, -1.947815100e-05, 2.015720940e-08, -7.376117610e-12, -9.179351730e+02, 6.830102380e-01, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.337279200e+00, -4.940247310e-05, 4.994567780e-07, -1.795663940e-10, 2.002553760e-14, -9.501589220e+02, -3.205023310e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
