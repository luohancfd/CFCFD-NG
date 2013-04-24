-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

NO = {}
NO.M = {
   value = 3.0006100e-02,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}
NO.atomic_constituents = {N=1,O=1}
NO.charge = 0
NO.gamma = {
   value = 1.3859e+00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/(Cp - R) from Chemkin-II coefficients'
}
NO.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 4.218476300e+00, -4.638976000e-03, 1.104102200e-05, -9.336135400e-09, 2.803577000e-12, 9.844623000e+03, 2.280846400e+00, }
   },
   { T_low  = 1000.0,
     T_high = 6000.0,
     coeffs = {0.000000000e+00, 0.000000000e+00, 3.260605600e+00, 1.191104300e-03, -4.291704800e-07, 6.945766900e-11, -4.033609900e-15, 9.920974600e+03, 6.369302700e+00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
