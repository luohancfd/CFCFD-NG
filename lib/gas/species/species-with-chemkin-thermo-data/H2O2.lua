-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

H2O2 = {}
H2O2.M = {
   value = 3.4014680e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
H2O2.gamma = {
   value = 1.2435e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
H2O2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.276112690e+00, -5.428224170e-04, 1.673357010e-05, -2.157708130e-08, 8.624543630e-12, -1.770258210e+04, 3.435050740e+00, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.165002850e+00, 4.908316940e-03, -1.901392250e-06, 3.711859860e-10, -2.879083050e-14, -1.786178770e+04, 2.916156620e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
