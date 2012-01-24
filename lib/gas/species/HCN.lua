-- Collater: Rowan J. Gollan
-- Date: 21-Jun-2009

HCN = {}
HCN.M = {
   value = 27.0253400e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
HCN.gamma = {
   value = 1.301,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
HCN.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = {  9.098286930e+04, -1.238657512e+03,  8.721307870e+00,
               -6.528242940e-03,  8.872700830e-06, -4.808886670e-09,
                9.317898500e-13,  2.098915450e+04, -2.746678076e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.236889278e+06, -4.446732410e+03,  9.738874850e+00,
               -5.855182640e-04,  1.072791440e-07, -1.013313244e-11,
                3.348247980e-16,  4.221513770e+04, -4.005774072e+01
    }
  },
  ref="from CEA2::thermo.inp"
}

HCN.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=300.0, T_high=1000.0, A=0.94863717e+00, B=-0.14891490e+03, C=0.15258721e+05, D=-0.72592817e+00},
      {T_low=1000.0, T_high=5000.0, A=0.57370725e+00, B=-0.85239973e+03, C=0.17953641e+06, D=0.24032031e+01},
      ref = 'from CEA2::trans.inp which cites Zeleznik & Svehla (1970), Svehla (1994)'
   }
}
HCN.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=300.0, T_high=1000.0, A=0.11749061e+01, B=-0.19100307e+03, C=0.15714065e+05, D=-0.13488014e+01},
      {T_low=1000.0, T_high=5000.0, A=0.50543688e+00, B=-0.13891056e+04, C=0.28003144e+06, D=0.42095130e+01},
      ref = 'from CEA2::trans.inp which cites Zeleznik & Svehla (1970), Svehla (1994)'
   }
}
