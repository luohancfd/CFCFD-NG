-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C2H = {}
C2H.M = {
   value = 25.029340e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H.gamma = {
   value = 1.2465e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 2.889657330e00, 1.340996110e-02, -2.847695010e-05, 2.947910450e-08, -1.093315110e-11, 6.683939320e04, 6.222964380e00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e00, 0.000000000e00, 3.167806520e00, 4.752219020e-03, -1.837870770e-06, 3.041902520e-10, -1.772327700e-14, 6.712106500e04, 6.635894750e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
