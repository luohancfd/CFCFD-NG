CI = {
  i = N,
  j = N,
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0101,   0.2216,  -1.8319,   8.3527 },
      Pi_Omega_22 = {  -0.0082,   0.1830,  -1.5724,   7.8962 },
      D           = {   0.0101,  -0.2216,   3.3319, -11.0152 },
    }
  }
}
