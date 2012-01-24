CI = {
  i = e_minus,
  j = O2,
  reference = 'Wright et al, AIAA Journal Vol. 43 No. 12 December 2005',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =  1000.0,
      T_high = 20000.0,
      Pi_Omega_11 = {   0.0175,  -0.5440,   5.5384, -16.6444 },
      Pi_Omega_22 = {   0.0097,  -0.3393,   3.8181, -11.9787 },
      D           = {  -0.0175,   0.5440,  -4.0384,  18.7092 },
    }
  }
}
