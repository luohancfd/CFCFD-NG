"""
least_squares.py: Fits a generalized polynomial basis to given data.

.. Author: Rowan J. Gollan
.. Version: 11-May-04
.. Version: 26-Jan-2012 PJ docs and example code

Example transcript::

    $python ~/e3bin/cfpylib/nm/least_squares.py
    Should be roughly y = x**2, i.e. coefficients [0,0,1].
    coeff= [0.06547619047616493, -0.092142857142826, 1.038095238095233]
"""

import math
from matrix import *

def least_squares(x, y, fList) :
    """ 
    Fit a generalized polynomial to the given data (x, y).

    The fitted model is of the form:
    y(x) = a[0]*fList[0](x) + a[1]*fList[1](x) + ... + a[m]*fList[m](x)

    :param x: list or array of data for independent variable
    :param y: list or array of data for dependent variable
    :param fList: list or array of basis functions to use

    :returns: a vector of coefficients a[]
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

if __name__ == '__main__':
    basis = [lambda x: 1.0, lambda x: x, lambda x: x*x]
    x_data = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    y_data = [0.1, 0.99, 3.9, 9.2, 16.6, 25.2, 37.0]
    coeff = least_squares(x_data, y_data, basis)
    print "Should be roughly y = x**2, i.e. coefficients [0,0,1]."
    print "coeff=", coeff
