-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

H2O = {}
H2O.M = {
   value = 18.01528e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
H2O.gamma = {
   value = 1.329,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
H2O.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=373.2, T_high=1073.2, A=0.50019557e+00, B=-0.69712796e+03, C=0.88163892e+05, D=0.30836508e+01},
      {T_low=1073.2, T_high=5000.0, A=0.58988538e+00, B=-0.53769814e+03, C=0.54263513e+05, D=0.23386375e+01},
      {T_low=5000.0, T_high=15000.0, A=0.64330087e+00, B=-0.95668913e+02, C=-0.37742283e+06, D=0.18125190e+01},
      ref = 'from CEA2::trans.inp which cites Sengers and Watson (1986)'
   }
}
H2O.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=373.2, T_high=1073.2, A=0.10966389e+01, B=-0.55513429e+03, C=0.10623408e+06, D=-0.24664550e+00},
      {T_low=1073.2, T_high=5000.0, A=0.39367933e+00, B=-0.22524226e+04, C=0.61217458e+06, D=0.58011317e+01},
      {T_low=5000.0, T_high=15000.0, A=-0.41858737e+00, B=-0.14096649e+05, C=0.19179190e+08, D=0.14345613e+02},
      ref = 'from CEA2::trans.inp which cites Sengers and Watson (1986)'
   }
}
H2O.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -3.947960830e+04,  5.755731020e+02,  9.317826530e-01,
	        7.222712860e-03, -7.342557370e-06,  4.955043490e-09,
	       -1.336933246e-12, -3.303974310e+04,  1.724205775e+01
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.034972096e+06, -2.412698562e+03,  4.646110780e+00,
                2.291998307e-03, -6.836830480e-07,  9.426468930e-11,
               -4.822380530e-15, -1.384286509e+04, -7.978148510e+00
    }
  },
  ref="from CEA2::thermo.inp"
}
H2O.T_c = {
   value = 647.14,
   units = 'K',
   description = 'critical temperature',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
H2O.p_c = {
   value = 220.64e+05,
   units = 'Pa',
   description = 'critical pressure',
   reference = 'Poling, B.E. et al. (2001). The Properties of Gases and Liquids. Section A, p.A.5'
}
