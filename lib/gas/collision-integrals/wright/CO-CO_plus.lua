CI = {
  i = 'CO',
  j = 'CO_plus',
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 15000.0,
      Pi_Omega_11 = {  -0.0025,   0.0975,  -1.2366,   9.1939 },
      Pi_Omega_22 = {  -0.0000,   0.0004,  -0.5035,   7.6373 },
      D           = {   0.0025,  -0.0975,   2.7366, -12.2029 },
    }
  }
}
