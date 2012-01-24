-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

NO2 = {}
NO2.M = {
   value = 46.0055e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
NO2.gamma = {
   value = 1.287,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
NO2.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=300.0, T_high=1000.0, A=0.57379100e+00, B=-0.12636034e+03, C=0.21566823e+04, D=0.22287492e+01},
      {T_low=1000.0, T_high=5000.0, A=0.64239645e+00, B=0.60012144e+00, C=-0.27020876e+05, D=0.16570566e+01},
      ref = 'from CEA2::trans.inp which cites Svehla (1966)'
   }
}
NO2.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=300.0, T_high=1000.0, A=0.48574998e+00, B=-0.50702110e+03, C=0.46605820e+05, D=0.36444556e+01},
      {T_low=1000.0, T_high=5000.0, A=0.97660465e+00, B=0.72760751e+03, C=-0.32527989e+06, D=-0.60899123e+00},
      ref = 'from CEA2::trans.inp which cites Bousheri et al (1987), Svehla (1966)'
   }
}
NO2.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -5.642038780e+04,  9.633085720e+02, -2.434510974e+00,
                1.927760886e-02, -1.874559328e-05,  9.145497730e-09,
               -1.777647635e-12, -1.547925037e+03,  4.067851210e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  7.213001570e+05, -3.832615200e+03,  1.113963285e+01,
               -2.238062246e-03,  6.547723430e-07, -7.611335900e-11,
                3.328361050e-15,  2.502497403e+04, -4.305130040e+01
    }
  },
  ref="from CEA2::thermo.inp"
}
