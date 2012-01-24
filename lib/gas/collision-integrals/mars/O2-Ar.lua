CI = {
  i = O2,
  j = Ar,
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 15000.0,
      Pi_Omega_11 = {  -0.0063,   0.1430,  -1.2832,   7.4047 },
      Pi_Omega_22 = {  -0.0089,   0.1990,  -1.6574,   8.2803 },
      D           = {   0.0063,  -0.1430,   2.7832, -10.5327 },
    }
  }
}
