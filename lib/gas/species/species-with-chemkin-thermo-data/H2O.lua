-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

H2O = {}
H2O.M = {
   value = 1.8015280e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
H2O.atomic_constituents = {H=2,O=1}
H2O.charge = 0
H2O.gamma = {
   value = 1.3289e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
H2O.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.198640560e+00, -2.036434100e-03, 6.520402110e-06, -5.487970620e-09, 1.771978170e-12, -3.029372670e+04, -8.490322080e-01, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.033992490e+00, 2.176918040e-03, -1.640725180e-07, -9.704198700e-11, 1.682009920e-14, -3.000429710e+04, 4.966770100e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
