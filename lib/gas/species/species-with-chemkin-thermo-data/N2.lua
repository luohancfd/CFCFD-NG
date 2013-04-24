-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

N2 = {}
N2.M = {
   value = 2.8013400e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
N2.atomic_constituents = {N=2}
N2.charge = 0
N2.gamma = {
   value = 1.4005e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
N2.CEA_coeffs = {
   { T_low  = 300.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.298677000e+00, 1.408240400e-03, -3.963222000e-06, 5.641515000e-09, -2.444854000e-12, -1.020899900e+03, 3.950372000e+00, }
   },
   { T_low  = 1000.0,
     T_high = 5000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 2.926640000e+00, 1.487976800e-03, -5.684760000e-07, 1.009703800e-10, -6.753351000e-15, -9.227977000e+02, 5.980528000e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
