-- Collater: Daniel F. Potter
-- Date: 15 July 2012

-- Ground state of atomic Carbon, 2s22p2 3P

C_3P = {}
C_3P.M = { 
   value = 12.0107000e-3,
   units = 'kg/mol',
   description = 'molecular mass',
   reference = 'molecular weight from CEA2'
}

-- Nonequilibrium data

C_3P.species_type = "monatomic"
C_3P.h_f = {
   value = 59670127.47,
   units = 'J/kg',
   description = 'Heat of formation',
   reference = 'from CEA2::thermo.inp'
}
C_3P.Z = {
   value = 0,
   units = 'ND',
   description = 'Charge number',
   reference = 'NA'
}
C_3P.electronic_levels = {
   n_levels = 1,
   ref = 'NIST ASD: http://physics.nist.gov/PhysRefData/ASD/index.html',
   -- ===========================================================
   --   No.      n       E(cm-1)     g     l     L     S     parity 
   -- ===========================================================
   ilev_0   =  { 2,       29.58,     9,   -1,    1,    1,    2 }
   -- ===========================================================
}

