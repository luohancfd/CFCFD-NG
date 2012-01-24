-- Collater: Brendan T. O'Flaherty
-- Date: 07 Aug 2009

C2H2 = {}
C2H2.M = {
   value = 26.037280e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}
C2H2.gamma = {
   value = 1.2321e00,
   units = 'non-dimensional',
   description = 'ratio of specific heats at 300.0K',
   reference = 'evaluated using Cp/R from Chemkin-II coefficients'
}
C2H2.CEA_coeffs = {
   { T_low  = 200.0,
     T_high = 1000.0,
     coeffs = {0.000000000e00, 0.000000000e00, 8.086810940e-01, 2.336156290e-02, -3.551718150e-05, 2.801524370e-08, -8.500729740e-12, 2.642898070e04, 1.393970510e01, }
   },
   { T_low  = 1000.0,
     T_high = 3500.0,
     coeffs = {0.000000000e00, 0.000000000e00, 4.147569640e00, 5.961666640e-03, -2.372948520e-06, 4.674121710e-10, -3.612352130e-14, 2.593599920e04, -1.230281210e00, }
   },
   ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
}
C2H2.T_c = {
   value = 308.30,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
C2H2.p_c = {
   value = 61.14e05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
