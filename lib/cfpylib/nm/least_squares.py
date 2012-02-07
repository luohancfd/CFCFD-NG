"""
least_squares.py: Fits a least-squares polynomial to given data.

Author: Rowan J. Gollan
Version: 11-May-04
"""

import math
from matrix import *

def least_squares(x, y, fList) :
    """ Returns a least-squares fit to the given data (x, y)
    as a polynomial of order, m.

    :param x: independent variable
    :param y: dependent variable
    :param m: order of least-squarse polynomial

    :returns: a vector of polynomial coefficients
    """
    if len(x) != len(y) :
        raise IndexError, 'Vectors of data (x, y) do not match in length'
            
    nr, nc = len(x), len(fList)
    phi = new_matrix(nr, nc)

    # Now fill in phi
    for i in range(nr) : # loops through the rows
        for j in range(len(fList)) : # loops through the columns
            phi[i][j] = fList[j](x[i])

    # Form the transpose of phi
    phiT = transpose(phi)

    # Find our matrix system to solve in the form
    # A x = b
    A = mmult(phiT, phi)

    # y is only a list, make it a vector
    yvec = [ y ]

    b = mmult(phiT, transpose(yvec))

    a = gauss_jordan(A, b, 0)
    # Now we tidy up the list format so that
    # the user just a gets a list of the coefficients

    a = transpose(a)
    a = a[0]

    return a
