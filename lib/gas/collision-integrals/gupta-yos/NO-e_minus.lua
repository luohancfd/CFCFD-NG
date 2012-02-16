CI = {
  i = ''NO'',
  j = ''e_minus'',
  reference = 'NASA Reference Publication 1232 by Gupta,R.N., Yos,J.M. et al 1990',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high =  8000.0,
      Pi_Omega_11 = {  -0.2202,    5.2265,  -40.5659,  104.7126, },
      Pi_Omega_22 = {  -0.2202,    5.2265,  -40.5659,  104.7126, },
      D           = {   0.2202,   -5.2261,   42.0630, -106.0937, },
    },
    {
      T_low  =  8000.0,
      T_high = 30000.0,
      Pi_Omega_11 = {  -0.2871,    8.3757,  -81.3787,  265.6292, },
      Pi_Omega_22 = {  -0.2871,    8.3757,  -81.3787,  265.6292, },
      D           = {   0.2871,   -8.3759,   82.8802, -267.0227, },
    },
  },
}
