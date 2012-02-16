CI = {
  i = 'CO2',
  j = 'N',
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0075,   0.1839,  -1.6674,   8.5548 },
      Pi_Omega_22 = {  -0.0085,   0.2066,  -1.8269,   8.9958 },
      D           = {   0.0075,  -0.1839,   3.1674, -11.4256 },
    }
  }
}
