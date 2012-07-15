-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 2nd excited state of nitrogen, B3Pig

N2_B = {}
N2_B.M = { 
   value = 28.0134e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'from CEA2::thermo.inp'
}

-- Nonequilibrium data

N2_B.species_type = "nonpolar diatomic"
N2_B.h_f = {
   value = 0.0,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
N2_B.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
N2_B.electronic_levels = {
   n_levels = 1,
   ref = 'Spradian07::diatom.dat',
   -- ===========================================================================================================================================================
   --   n       Te         re       g   dzero      we         wexe      weye        weze        be        alphae      de          betae       spn-orb     l   s  
   -- ===========================================================================================================================================================
   ilev_0  = {  59619.09,  1.2126,  6,  38660.00,  1734.025,  14.4120, -3.300E-03, -7.900E-04,  1.63772,  1.793E-02,  5.900E-06,  0.000E+00,  4.224E+01,  1,  3 }
   -- ===========================================================================================================================================================
}
   
