CI = {
  i = 'e_minus',
  j = 'N2',
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   500.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0290,   0.5808,  -3.4068,   7.3027 },
      Pi_Omega_22 = {  -0.0303,   0.5681,  -2.9363,   5.0715 },
      D           = {   0.0290,  -0.5808,   4.9068,  -5.2378 },
    }
  }
}
