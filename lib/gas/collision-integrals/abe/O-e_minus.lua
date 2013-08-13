CI = {
  i = 'O',
  j = 'e_minus',
  reference = 'NASA Reference Publication 1232 by Gupta,R.N., Yos,J.M. et al 1990',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high =  9000.0,
      Pi_Omega_11 = {   0.0164,   -0.2431,    1.1231,   -1.5561, },
      Pi_Omega_22 = {   0.0164,   -0.2431,    1.1231,   -1.5561, },
      D           = {   0.0581,   -1.5975,   15.4508,  -40.7370, },
    },
    {
      T_low  =  9000.0,
      T_high = 30000.0,
      Pi_Omega_11 = {  -0.2027,    5.6428,  -51.5646,  155.6091, },
      Pi_Omega_22 = {  -0.2027,    5.6428,  -51.5646,  155.6091, },
      D           = {   0.0581,   -1.5975,   15.4508,  -40.7370, },
    },
  },
}
