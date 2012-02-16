CI = {
  i = 'N2_plus',
  j = 'O2',
  reference = 'Wright et al, JTHT Vol. 19 No. 1 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 12000.0,
      Pi_Omega_11 = {   0.0088,  -0.1140,  -0.2169,   7.9054 },
      Pi_Omega_22 = {   0.0101,  -0.1415,   0.0087,   7.3237 },
      D           = {  -0.0088,   0.1140,   1.7169, -10.9465 },
    }
  }
}
