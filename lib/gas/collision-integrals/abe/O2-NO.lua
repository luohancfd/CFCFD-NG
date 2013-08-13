CI = {
  i = 'O2',
  j = 'NO',
  reference = 'NASA Reference Publication 1232 by Gupta,R.N., Yos,J.M. et al 1990',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high = 30000.0,
      Pi_Omega_11 = {   0.0000,   -0.0438,    0.5352,    1.7252, },
      Pi_Omega_22 = {   0.0000,   -0.0522,    0.7045,    1.0738, },
      D           = {   0.0000,    0.0438,    0.9647,   -8.2380, },
    },
  },
}
