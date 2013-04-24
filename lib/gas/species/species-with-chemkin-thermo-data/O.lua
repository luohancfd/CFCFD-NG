-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

O = {}
O.M = {
   value = 1.5999400e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
O.atomic_constituents = {O=1}
O.charge = 0
O.gamma = {
   value = 1.6120e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
O.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.168267100e+00, -3.279318840e-03, 6.643063960e-06, -6.128066240e-09, 2.112659710e-12, 2.912225920e+04, 2.051933460e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.569420780e+00, -8.597411370e-05, 4.194845890e-08, -1.001777990e-11, 1.228336910e-15, 2.921757910e+04, 4.784338640e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
