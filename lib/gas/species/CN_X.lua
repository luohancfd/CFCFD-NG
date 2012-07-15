-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- Ground state of cyanogen, X2Sigma+

CN_X = {}
CN_X.M = {
   value = 26.0174000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}

-- Nonequilibrium data

CN_X.species_type = "polar diatomic"
CN_X.h_f = {
   value = 16861160.30,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CN_X.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CN_X.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.1718,  2,  61631.00,  2068.640,  13.0988, -1.130E-02,  3.200E-04,  1.89978,  1.736E-02,  6.400E-06,  1.200E-08,  0.000E+00,  0,  2 }
   -- ===========================================================================================================================================================
}
