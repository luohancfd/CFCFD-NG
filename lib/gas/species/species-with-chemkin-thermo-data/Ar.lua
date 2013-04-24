-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

Ar = {}
Ar.M = {
   value = 3.9948000e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
Ar.atomic_constituents = {Ar=1}
Ar.charge = 0
Ar.gamma = {
   value = 1.6667e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
Ar.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, -7.453750000e+02, 4.366000000e+00, }
   },
   { T_low  = 1000.0,
     T_high = 5000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, 0.000000000e+00, -7.453750000e+02, 4.366000000e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
