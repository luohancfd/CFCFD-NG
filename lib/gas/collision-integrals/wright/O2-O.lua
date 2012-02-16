CI = {
  i = 'O2',
  j = 'O',
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 10000.0,
      Pi_Omega_11 = {  -0.0055,   0.1174,  -1.0770,   6.6923 },
      Pi_Omega_22 = {  -0.0049,   0.1020,  -0.9228,   6.3111 },
      D           = {   0.0055,  -0.1174,   2.5770,  -9.5651 },
    }
  }
}
