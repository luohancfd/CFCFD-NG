CI = {
  i = CO2,
  j = O,
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   500.0,
      T_high = 20000.0,
      Pi_Omega_11 = {   0.0002,  -0.0134,  -0.0206,   3.9632 },
      Pi_Omega_22 = {  -0.0009,   0.0120,  -0.2145,   4.5663 },
      D           = {  -0.0002,   0.0134,   1.5206,  -6.8836 },
    }
  }
}
