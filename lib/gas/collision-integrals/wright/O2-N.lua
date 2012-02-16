CI = {
  i = 'O2',
  j = 'N',
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   500.0,
      T_high = 15000.0,
      Pi_Omega_11 = {  -0.0014,   0.0188,  -0.2741,   4.4792 },
      Pi_Omega_22 = {  -0.0010,   0.0104,  -0.2036,   4.4182 },
      D           = {   0.0014,  -0.0188,   1.7741,  -7.3067 },
    }
  }
}
