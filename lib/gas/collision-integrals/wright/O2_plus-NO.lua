CI = {
  i = O2_plus,
  j = NO,
  reference = 'Wright et al, JTHT Vol. 19 No. 1 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 12000.0,
      Pi_Omega_11 = {   0.0291,  -0.5849,   3.2164,  -0.1550 },
      Pi_Omega_22 = {   0.0398,  -0.8408,   5.2681,  -5.5231 },
      D           = {  -0.0291,   0.5849,  -1.7164,  -2.9042 },
    }
  }
}
