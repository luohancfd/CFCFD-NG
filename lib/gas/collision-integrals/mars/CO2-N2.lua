CI = {
  i = 'CO2',
  j = 'N2',
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0077,   0.1751,  -1.5137,   8.1593 },
      Pi_Omega_22 = {  -0.0119,   0.2770,  -2.3225,  10.3644 },
      D           = {   0.0077,  -0.1751,   3.0137, -11.2686 },
    }
  }
}
