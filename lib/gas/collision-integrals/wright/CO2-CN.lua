CI = {
  i = CO2,
  j = CN,
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   500.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0008,   0.0091,  -0.1967,   4.7231 },
      Pi_Omega_22 = {  -0.0008,   0.0102,  -0.1988,   4.8355 },
      D           = {   0.0008,  -0.0091,   1.6967,  -7.8095 },
    }
  }
}
