-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 2nd excited state of cyanogen, B2Sigma+

CN_B = {}
CN_B.M = {
   value = 26.0174000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}

-- Nonequilibrium data

CN_B.species_type = "polar diatomic"
CN_B.h_f = {
   value = 16861160.30,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CN_B.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CN_B.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  25752.80,  1.1511,  2,  55130.00,  2161.420,  18.1200, -4.309E-01,  0.000E+00,  1.96882,  2.005E-02,  6.600E-06,  0.000E+00,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
