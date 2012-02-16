CI = {
  i = 'O2',
  j = 'O2',
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 15000.0,
      Pi_Omega_11 = {  -0.0023,   0.0516,  -0.5785,   5.6041 },
      Pi_Omega_22 = {  -0.0089,   0.2066,  -1.7522,   8.6099 },
      D           = {   0.0023,  -0.0516,   2.0785,  -8.6796 },
    }
  }
}
