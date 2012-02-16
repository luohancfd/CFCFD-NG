CI = {
  i = 'CO2',
  j = 'O2',
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0068,   0.1637,  -1.4935,   8.2050 },
      Pi_Omega_22 = {  -0.0110,   0.2678,  -2.3236,  10.4795 },
      D           = {   0.0068,  -0.1637,   2.9935, -11.3539 },
    }
  }
}
