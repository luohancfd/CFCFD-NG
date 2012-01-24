CI = {
  i = CO2,
  j = CO,
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0078,   0.1780,  -1.5361,   8.2160 },
      Pi_Omega_22 = {  -0.0119,   0.2770,  -2.3226,  10.3647 },
      D           = {   0.0078,  -0.1780,   3.0361, -11.3252 },
    }
  }
}
