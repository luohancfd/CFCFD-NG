CI = {
  i = 'e_minus',
  j = 'O',
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high = 20000.0,
      Pi_Omega_11 = {   0.0085,  -0.2361,   2.3343,  -6.8487 },
      Pi_Omega_22 = {  -0.0193,   0.4372,  -2.9464,   6.7816 },
      D           = {  -0.0085,   0.2361,  -0.8343,   8.9136 },
    }
  }
}
