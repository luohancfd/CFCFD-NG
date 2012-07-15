-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- 3rd excited state of atomic Carbon, 2s2p3 5S¡

C_5So = {}
C_5So.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}

-- Nonequilibrium data

C_5So.species_type = "monatomic"
C_5So.h_f = {
   value = 59670127.47,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C_5So.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C_5So.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n       E(cm-1)     g     l     L     S     parity 
   -- ===========================================================
   ilev_0   =  { 2,    33735.20,     5,   -1,    0,    2,    1 }
   -- ===========================================================
}

