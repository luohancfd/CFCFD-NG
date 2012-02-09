import numpy as np

# Some constants

_C = 2.99792458e8
_EMASS = 511.003414e3


def calcB0(B_Rem, lambda_0, gap):
    """Calculate B0 for an undulator of wavelength lambda at gap gap"""
    return B_rem / np.sinh(np.pi * gap / lambda_0)

def calcK(lambda_0, B):
    """Calculate K for an undulator"""
    return _C * lambda_0 * B / (_EMASS * 2.0 * np.pi)


