CI = {
  i = 'O2_plus',
  j = 'O2',
  reference = 'Wright et al, JTHT Vol. 19 No. 1 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 12000.0,
      Pi_Omega_11 = {  -0.0107,   0.2670,  -2.3710,  11.8459 },
      Pi_Omega_22 = {   0.0410,  -0.8580,   5.3038,  -5.3605 },
      D           = {   0.0107,  -0.2670,   3.8710, -14.9214 },
    }
  }
}
