-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 5th excited state of carbon monoxide, A1Pi

CO_A = {}
CO_A.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}

-- Nonequilibrium data

CO_A.species_type = "polar diatomic"
CO_A.h_f = {
   value = -3946262.10,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CO_A.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CO_A.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  65075.77,  1.2353,  2,  24740.00,  1518.240,  19.4000,  7.660E-01,  0.000E+00,  1.61150,  2.325E-02,  7.330E-06,  1.000E-07,  0.000E+00,  1,  1 }
   -- ===========================================================================================================================================================
}
