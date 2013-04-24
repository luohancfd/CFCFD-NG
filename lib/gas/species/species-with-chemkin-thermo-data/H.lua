-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

H = {}
H.M = {
   value = 1.0079400e-03,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
H.atomic_constituents = {H=1}
H.charge = 0
H.gamma = {
   value = 1.6667e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
H.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.500000000e+00, 7.053328190e-13, -1.995919640e-15, 2.300816320e-18, -9.277323320e-22, 2.547365990e+04, -4.466828530e-01, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.500000010e+00, -2.308429730e-11, 1.615619480e-14, -4.735152350e-18, 4.981973570e-22, 2.547365990e+04, -4.466829140e-01, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
