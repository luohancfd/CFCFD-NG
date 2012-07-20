-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of cyanogen, A2PiI

CN_A = {}
CN_A.M = {
   value = 26.0174000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'CEA2::thermo.inp'
}

-- Nonequilibrium data

CN_A.species_type = "polar diatomic"
CN_A.h_f = {
   value = 16861160.30,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CN_A.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CN_A.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {   9245.28,  1.2333,  4,  52580.00,  1812.560,  12.6090, -1.180E-02,  0.000E+00,  1.71510,  1.708E-02,  5.930E-06,  4.200E-08, -5.264E+01,  1,  2 }
   -- ===========================================================================================================================================================
}
