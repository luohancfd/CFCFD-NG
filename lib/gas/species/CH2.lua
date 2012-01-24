-- Collater: Rowan J. Gollan
-- Date: 20-Jun-2009

CH2 = {}
CH2.M = {
   value = 14.0265800e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
CH2.gamma = {
   value = 1.311,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
CH2.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  3.218921730e+04, -2.877601815e+02,  4.203583820e+00,
		3.455405960e-03, -6.746193340e-06,  7.654571640e-09,
	       -2.870328419e-12,  4.733624710e+04, -2.143628603e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  2.550418031e+06, -7.971625390e+03,  1.228924487e+01,
	       -1.699122922e-03,  2.991728605e-07, -2.767007492e-11,
                1.051341740e-15,  9.642216890e+04, -6.094739910e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
CH2.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.57643622e+00, B=-0.93704079e+02, C=0.86992395e+03, D=0.17333347e+01},
      {T_low=1000.0, T_high=5000.0, A=0.66400044e+00, B=0.10860843e+02, C=-0.76307841e+04, D=0.10323984e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}
CH2.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=200.0, T_high=1000.0, A=0.10238177e+01, B=-0.31092375e+03, C=0.32944309e+05, D=0.67787437e+00},
      {T_low=1000.0, T_high=5000.0, A=0.77485028e+00, B=-0.40089627e+03, C=-0.46551082e+05, D=0.25671481e+01},
      ref = 'ref="cea2 data for CH4"'
   }
}
