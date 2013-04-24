-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

HO2 = {}
HO2.M = {
   value = 3.3006740e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
HO2.atomic_constituents = {H=1,O=2}
HO2.charge = 0
HO2.gamma = {
   value = 1.3124e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
HO2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.301798010e+00, -4.749120510e-03, 2.115828910e-05, -2.427638940e-08, 9.292251240e-12, 2.948080400e+02, 3.716662450e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.017210900e+00, 2.239820130e-03, -6.336581500e-07, 1.142463700e-10, -1.079085350e-14, 1.118567130e+02, 3.785102150e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
