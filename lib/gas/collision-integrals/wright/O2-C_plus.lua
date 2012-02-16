CI = {
  i = 'O2',
  j = 'C_plus',
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {   0.0065,  -0.0774,  -0.4276,   8.2411 },
      Pi_Omega_22 = {   0.0133,  -0.2392,   0.8487,   5.0074 },
      D           = {  -0.0065,   0.0774,   1.9276, -11.0139 },
    }
  }
}
