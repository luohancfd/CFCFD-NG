-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

O2 = {}
O2.M = {
   value = 3.1998800e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
O2.gamma = {
   value = 1.3945e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
O2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.782456360e+00, -2.996734160e-03, 9.847302010e-06, -9.681295090e-09, 3.243728370e-12, -1.063943560e+03, 3.657675730e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.282537840e+00, 1.483087540e-03, -7.579666690e-07, 2.094705550e-10, -2.167177940e-14, -1.088457720e+03, 5.453231290e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
