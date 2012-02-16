CI = {
  i = 'O2',
  j = 'CO',
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0031,   0.0602,  -0.5824,   5.5388 },
      Pi_Omega_22 = {  -0.0079,   0.1733,  -1.4406,   7.7641 },
      D           = {   0.0031,  -0.0602,   2.0824,  -8.5799 },
    }
  }
}
