CI = {
  i = e_minus,
  j = CO2,
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high = 20000.0,
      Pi_Omega_11 = {   0.1041,  -2.3155,  16.0679, -31.4029 },
      Pi_Omega_22 = {   0.0580,  -1.1115,   5.8355,  -3.3769 },
      D           = {  -0.1041,   2.3155, -14.5679,  33.4678 },
    }
  }
}
