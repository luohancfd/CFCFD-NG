CI = {
  i = O2,
  j = Ar_plus,
  reference = 'Levin et al, JTHT Vol. 18 No. 1 2004',
  model = 'GuptaYos curve fits',
  parameters = {
    {
      T_low  =   300.0,
      T_high = 15000.0,
      Pi_Omega_11 = {   0.0080,  -0.1040,  -0.2728,   7.9528 },
      Pi_Omega_22 = {   0.0165,  -0.3037,   1.2944,   3.9660 },
      D           = {  -0.0080,   0.1040,   1.7728, -11.0808 },
    }
  }
}
