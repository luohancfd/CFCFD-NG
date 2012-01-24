CI = {
  i = N,
  j = Ar,
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   500.0,
      T_high = 15000.0,
      Pi_Omega_11 = {  -0.0011,   0.0100,  -0.2153,   4.4329 },
      Pi_Omega_22 = {  -0.0011,   0.0120,  -0.2241,   4.5798 },
      D           = {   0.0011,  -0.0100,   1.7153,  -7.2916 },
    }
  }
}
