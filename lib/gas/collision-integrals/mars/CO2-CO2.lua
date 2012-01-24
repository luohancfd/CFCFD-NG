CI = {
  i = CO2,
  j = CO2,
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0115,   0.2886,  -2.5580,  11.3378 },
      Pi_Omega_22 = {  -0.0148,   0.3731,  -3.2655,  13.3972 },
      D           = {   0.0115,  -0.2886,   4.0580, -14.5727 },
    }
  }
}
