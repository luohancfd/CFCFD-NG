-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 1st excited state of carbon monoxide, a 3Pir

CO_a3 = {}
CO_a3.M = {
   value = 28.010100e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'Periodic table'
}

-- Nonequilibrium data

CO_a3.species_type = "polar diatomic"
CO_a3.h_f = {
   value = -3946262.10,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
CO_a3.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
CO_a3.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  48686.70,  1.2057,  6,  41020.00,  1743.410,  14.3600, -4.500E-02,  0.000E+00,  1.69124,  1.904E-02,  6.360E-06,  4.000E-08,  4.153E+01,  1,  3 }
   -- ===========================================================================================================================================================
}
