CI = {
  i = 'O2',
  j = 'NO',
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 15000.0,
      Pi_Omega_11 = {  -0.0028,   0.0583,  -0.6053,   5.6543 },
      Pi_Omega_22 = {  -0.0099,   0.2246,  -1.8646,   8.8696 },
      D           = {   0.0028,  -0.0583,   2.1053,  -8.7135 },
    }
  }
}
