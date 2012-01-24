-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

N = {}
N.M = {
   value = 1.4006700e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
N.gamma = {
   value = 1.6667e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
N.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 5.610463700e+04, 4.193908700e+00, }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.415942900e+00, 1.748906500e-04, -1.190236900e-07, 3.022624500e-11, -2.036098200e-15, 5.613377300e+04, 4.649609600e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
