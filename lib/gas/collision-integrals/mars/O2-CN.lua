CI = {
  i = O2,
  j = CN,
  reference = 'Wright et al, AIAA Journal Vol. 45 No. 1 January 2007',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   500.0,
      T_high = 20000.0,
      Pi_Omega_11 = {  -0.0017,   0.0304,  -0.3720,   4.9929 },
      Pi_Omega_22 = {  -0.0015,   0.0252,  -0.3253,   4.9855 },
      D           = {   0.0017,  -0.0304,   1.8720,  -8.0140 },
    }
  }
}
