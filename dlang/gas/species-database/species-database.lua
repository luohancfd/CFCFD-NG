-- Curator: Rowan J. Gollan
-- Date: 08-Mar-2015
--
-- History:
--   08-Mar-2015 -- first cut.
--               -- Experiment to see if this form has longevity.
--

db = {}

db.N2 = {}
db.N2.atomic_constituents = {N=2}
db.N2.charge = 0
db.N2.M = { 
   value = 28.0134e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}
db.N2.gamma = { 
   value = 1.4,
   units = 'non-dimensional',
   description = '(ideal) ratio of specific heats at room temperature',
   reference = 'diatomic molecule at low temperatures, gamma = 7/5'
}
db.N2.entropyRefValues = {
   s1 = 0.0,
   T1 = 298.15,
   p1 = 101.325e3
}
db.N2.sutherlandVisc = {
   mu_ref = 1.663e-05,
   T_ref = 273.0,
   S = 107.0,
   ref = "Table 1-2, White (2006)"
}
db.N2.sutherlandThermCond = { 
   k_ref = 0.0242,
   T_ref = 273.0,
   S = 150.0,
   ref = "Table 1-3, White (2006)"
}
db.N2.ceaThermoCoeffs = {
   nsegments = 3,
   segment0 = {
      T_lower  = 200.0,
      T_upper = 1000.0,
      coeffs = { 2.210371497e+04,
		-3.818461820e+02,
		 6.082738360e+00,
		-8.530914410e-03,
		 1.384646189e-05,
	        -9.625793620e-09,
		 2.519705809e-12,
		 7.108460860e+02,
		-1.076003744e+01
      }
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 6000.0,
      coeffs = { 5.877124060e+05,
                -2.239249073e+03,
                 6.066949220e+00,
		-6.139685500e-04,
                 1.491806679e-07,
                -1.923105485e-11,
                 1.061954386e-15,
                 1.283210415e+04,
                -1.586640027e+01
      }
   },
   segment2 = {
      T_lower  = 6000.0,
      T_upper = 20000.0,
      coeffs = { 8.310139160e+08,
                -6.420733540e+05,
                 2.020264635e+02,
		-3.065092046e-02,
		 2.486903333e-06,
		-9.705954110e-11,
		 1.437538881e-15,
		 4.938707040e+06,
		-1.672099740e+03
      }
   },
   ref="from CEA2::thermo.inp"
}
db.N2.ceaViscosity = {
   nsegments = 3,
   segment0 = {
      A = 0.62526577,
      C = -1640.7983,
      B = -31.779652,
      T_upper = 1000,
      T_lower = 200,
      D = 1.7454992,
    },
   segment1 = {
      A = 0.87395209,
      C = -173948.09,
      B = 561.52222,
      T_upper = 5000,
      T_lower = 1000,
      D = -0.39335958,
    },
   segment2 = {
      A = 0.88503551,
      C = -731290.61,
      B = 909.02171,
      T_upper = 15000,
      T_lower = 5000,
      D = -0.53503838,
    },
}
db.N2.ceaThermCond = {
   nsegments = 3,
   segment0 = {
      A = 0.85439436,
      C = -12347.848,
      B = 105.73224,
      T_upper = 1000,
      T_lower = 200,
      D = 0.47793128,
    },
    segment1 =  {
      A = 0.88407146,
      C = -11429.64,
      B = 133.57293,
      T_upper = 5000,
      T_lower = 1000,
      D = 0.24417019,
    },
    segment2 = {
      A = 2.4176185,
      C = 3105580.2,
      B = 8047.7749,
      T_upper = 15000,
      T_lower = 5000,
      D = -14.517761,
    },
}
