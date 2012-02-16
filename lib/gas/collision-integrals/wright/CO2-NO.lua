CI = {
  i = 'CO2',
  j = 'NO',
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0072,   0.1707,  -1.5268,   8.2808 },
      Pi_Omega_22 = {  -0.0122,   0.2923,  -2.4941,  10.9021 },
      D           = {   0.0072,  -0.1707,   3.0268, -11.4108 },
    }
  }
}
