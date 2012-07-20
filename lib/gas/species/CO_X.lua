-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- Ground state of carbon monoxide, X 1Sigma+

CO_X = {}
CO_X.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}

-- Nonequilibrium data

CO_X.species_type = "polar diatomic"
CO_X.h_f = {
   value = -3946262.10,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CO_X.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CO_X.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {      0.00,  1.1283,  1,  89490.00,  2169.814,  13.2883,  1.051E-02,  0.000E+00,  1.93128,  1.750E-02,  6.121E-06, -1.153E-09,  0.000E+00,  0,  1 }
   -- ===========================================================================================================================================================
}
