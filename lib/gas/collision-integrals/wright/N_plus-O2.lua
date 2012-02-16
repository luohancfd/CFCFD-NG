CI = {
  i = 'N_plus',
  j = 'O2',
  reference = 'Wright et al, JTHT Vol. 19 No. 1 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 12000.0,
      Pi_Omega_11 = {   0.0199,  -0.3613,   1.4902,   4.1244 },
      Pi_Omega_22 = {   0.0286,  -0.5693,   3.1572,  -0.2295 },
      D           = {  -0.0199,   0.3613,   0.0098,  -6.9519 },
    }
  }
}
