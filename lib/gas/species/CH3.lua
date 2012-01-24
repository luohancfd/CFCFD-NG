-- Collater: Rowan J. Gollan
-- Date: 17-Apr-2009

CH3 = {}
CH3.M = {
   value = 15.0345200e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CH3.gamma = {
   value = 1.276,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CH3.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -2.876188806e+04,  5.093268660e+02,  2.002143949e-01,
	        1.363605829e-02, -1.433989346e-05,  1.013556725e-08,
               -3.027331936e-12,  1.408271825e+04,  2.022772791e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  2.760802663e+06, -9.336531170e+03,  1.487729606e+01,
               -1.439429774e-03,  2.444477951e-07, -2.224555778e-11,
                8.395065760e-16,  7.481809480e+04, -7.919682400e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
-- CH3.CEA_coeffs = {
--    { T_low  = 200.0,
--      T_high = 1000.0,
--      coeffs = {-2.87618881e+04, 5.09326866e+02, 2.00214395e-01, 1.36360583e-02, -1.43398935e-05, 1.01355673e-08, -3.02733194e-12, 1.40827182e+04, 2.02277279e+01}
--    },
--    { T_low  = 1000.0,
--      T_high = 6000.0,
--      coeffs = {2.76080266e+06, -9.33653117e+03, 1.48772961e+01, -1.43942977e-03, 2.44447795e-07, -2.22455578e-11, 8.39506576e-16, 7.48180948e+04, -7.91968240e+01}
--    },
--    ref='The Chemkin Thermodynamic Data Base, Kee R.J. et al (1993)'
-- }
CH3.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.57643622e+00, B=-0.93704079e+02, C=0.86992395e+03, D=0.17333347e+01},
      {T_low=1000.0, T_high=5000.0, A=0.66400044e+00, B=0.10860843e+02, C=-0.76307841e+04, D=0.10323984e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}
CH3.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.10238177e+01, B=-0.31092375e+03, C=0.32944309e+05, D=0.67787437e+00},
      {T_low=1000.0, T_high=5000.0, A=0.77485028e+00, B=-0.40089627e+03, C=-0.46551082e+05, D=0.25671481e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}
