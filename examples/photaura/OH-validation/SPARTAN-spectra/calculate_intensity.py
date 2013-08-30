from cfpylib.util.YvX import YvX
from math import exp
from radpy import *

dL = 0.1	# path length in metres
HWHM = 2.5	# HWHM of apparatus function in Ang


JvE = YvX("nu_IE_IA_69_no_AF.txt",0,1)
KvE = YvX("nu_IE_IA_69_no_AF.txt",0,2)

S = SpectralIntensity()

for i,eta in enumerate(JvE.x_array):
    kappa_eta = KvE.y_array[i]
    j_eta = JvE.y_array[i]
    if kappa_eta==0.0:
        I_eta = j_eta * dL
    else:
        I_eta = j_eta / kappa_eta * ( 1.0 - exp( - dL * kappa_eta ) )
    nu = eta * RC_c
    S.nu.push_back(nu)
    S.I_nu.push_back(I_eta * eta / nu)

S.I_int.resize( S.nu.size() )

# Apparatus function
A = Voigt(0.0, HWHM, 1)
S.apply_apparatus_function(A)

S.write_to_file("intensity_spectra.txt")
