CI = {
  i = ''O2'',
  j = ''e_minus'',
  reference = 'NASA Reference Publication 1232 by Gupta,R.N., Yos,J.M. et al 1990',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high =  9000.0,
      Pi_Omega_11 = {   0.0241,   -0.3467,    1.3887,   -0.0110, },
      Pi_Omega_22 = {   0.0241,   -0.3467,    1.3887,   -0.0110, },
      D           = {  -0.0241,    0.3464,    0.1136,   -1.3848, },
    },
    {
      T_low  =  9000.0,
      T_high = 30000.0,
      Pi_Omega_11 = {   0.0025,   -0.0742,    0.7235,   -0.2116, },
      Pi_Omega_22 = {   0.0025,   -0.0742,    0.7235,   -0.2116, },
      D           = {  -0.0029,    0.0856,    0.6655,   -0.8205, },
    },
  },
}
