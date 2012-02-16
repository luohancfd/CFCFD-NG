CI = {
  i = 'O2_plus',
  j = 'N2',
  reference = 'Wright et al, JTHT Vol. 19 No. 1 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 12000.0,
      Pi_Omega_11 = {   0.0273,  -0.5379,   2.8060,   1.0689 },
      Pi_Omega_22 = {   0.0399,  -0.8424,   5.2572,  -5.3935 },
      D           = {  -0.0273,   0.5379,  -1.3060,  -4.1100 },
    }
  }
}
