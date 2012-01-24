-- Collater: Rowan J. Gollan
-- Date: 30-Mar-2009

OH = {}
OH.M = {
   value = 17.00734e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}
OH.gamma = {
   value = 1.386,
   units = 'non-dimensional',
   description = 'ratio of specific heats at room temperature (= Cp/(Cp - R))',
   reference = 'using Cp evaluated from CEA2 coefficients at T=300.0 K'
}
OH.viscosity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.59711536e+00, B=-0.46100678e+03, C=0.37606286e+05, D=0.24041761e+01},
      {T_low=5000.0, T_high=15000.0, A=0.64287721e+00, B=-0.18173747e+03, C=-0.88543767e+05, D=0.19636057e+01},
      ref = 'from CEA2::trans.inp which cites Svehla (1994)'
   }
}
OH.thermal_conductivity = {
   model = "CEA",
   parameters = {
      {T_low=1000.0, T_high=5000.0, A=0.68627561e+00, B=-0.74033274e+03, C=0.27559033e+05, D=0.28308741e+01},
      {T_low=5000.0, T_high=15000.0, A=-0.47918112e+00, B=-0.93769908e+04, C=0.70509952e+07, D=0.14203688e+02},
      ref = 'from CEA2::trans.inp which cites Svehla (1994)'
   }
}
OH.CEA_coeffs = {
  { T_low  = 200.0,
    T_high = 1000.0,
    coeffs = { -1.998858990e+03,  9.300136160e+01,  3.050854229e+00,
	        1.529529288e-03, -3.157890998e-06,  3.315446180e-09,
	       -1.138762683e-12,  2.991214235e+03,  4.674110790e+00
    }
  },
  { T_low  = 1000.0,
    T_high = 6000.0,
    coeffs = {  1.017393379e+06, -2.509957276e+03,  5.116547860e+00,
		1.305299930e-04, -8.284322260e-08,  2.006475941e-11,
               -1.556993656e-15,  2.019640206e+04, -1.101282337e+01
    }
  },
  { T_low = 6000.0,
    T_high = 20000.0,
    coeffs = {  2.847234193e+08, -1.859532612e+05,  5.008240900e+01,
               -5.142374980e-03,  2.875536589e-07, -8.228817960e-12,
                9.567229020e-17,  1.468393908e+06, -4.023555580e+02
    }
 },
  ref="from CEA2::thermo.inp"
}
