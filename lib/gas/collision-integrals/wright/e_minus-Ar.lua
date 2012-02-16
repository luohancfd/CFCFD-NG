CI = {
  i = 'e_minus',
  j = 'Ar',
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0175,   0.4452,  -2.7322,   2.2933 },
      Pi_Omega_22 = {  -0.1061,   2.5811, -19.6302,  46.6876 },
      D           = {   0.0175,  -0.4452,   4.2322,  -0.2285 },
    }
  }
}
