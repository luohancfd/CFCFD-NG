# file: ode_lists.py
"""Integrate a set of first-order ODEs using lists for storage.

This module provides a small number of ODE integration schemes
together with a coordination function ode_integrate().
The schemes will integrate a system of (nonstiff) ODEs
over a given range of the independent variable and with
a specified step.
The user is expected to supply a function that accepts the
current value of the independent variable and system state
and which provides the derivatives of the system as a list (or tuple).
For further details, see the doc string of ode_integrate() and study
the sample functions near the end of this file.

Running this module as a Python script gives me the following
transcript:
    [peterj@m4425 python_code]$ python ode.py
    Start sample integrations...
    (1) Constant derivatives with Euler method:
    x1= 10.0
    y1= [10.0, -20.0, 30.0]
    err1= [0.0, 0.0, 0.0]
    (2) Second-order linear ODE with Euler method:
    x1= 6.28318530718
    y1= [-0.00021600959761125483, 1.0319017729300233]
    err1= [0.020317859582370514, 0.020305117589602524]
    (3) Second-order linear ODE with RKF45 method:
    x1= 6.28318530718
    y1= [8.4505075463245372e-14, 1.000000000000574]
    err1= [5.1242191950857039e-11, 5.1282185272544397e-11]
    Done.

P. Jacobs
School of Engineering, UQ
05-Oct-03
"""

import math
import copy


def ode_integrate(x0, xlast, dx, dydx, n, y0, method="euler"):
    """Steps the set of ODEs until x reaches xlast.

    This function coordinates the work of integrating a system
    of first-order differential equations of the form
        y' = f(x, y)
        y(x=x0) = y0
    The actual work is done by one of a set of more specialised
    stepping functions that appear below.

    Input:
    x0    is the starting value of the independent variable
    xlast is the desired finishing value for x
    dx    is the requested step size
    dydx  is a callable function that returns the derivative of y wrt x
          The signature of this function is f(x, y, n) where
          x is a float value, y is a list (or array) or float values
          and n is an integer specifying the number of equations.
    n     is the number of dependent variables (in y)
    y0    is a list (or array) of starting values for the dependent variables
          It is assumed that the y-elements are indexed 0..n-1
    method is a string specifying which stepping method to use.

    Output:
    Final values of x, y, and error estimates for y values
    at x=xlast are returned as a tuple.    
    """
    assert callable(dydx)
    assert n <= len(y0)
    assert x0 <= xlast
    
    x = x0
    y = copy.copy(y0)  # Set up a new list so we don't mangle y0 itself.
    err_sum = [0.0] * n
    
    while x < xlast:
        # Decide if we need to take a full step
        if x + dx <= xlast:
            h = dx
        else:
            h = xlast - x
        # The real work is delegated.    
        if method == "euler":
            (x, y, err) = euler_step(x, h, dydx, n, y)
        elif method == "rkf45":
            (x, y, err) = rkf45_step(x, h, dydx, n, y)
        else:
            raise ValueError, "Requested stepping method is not available."
        for i in range(n): err_sum[i] += err[i]

    # Return the point reached, the set of y-values and an error estimate.
    return (x, y, err_sum)


def euler_step(x0, h, dydx, n, y0):
    "Steps the set of ODEs by the Euler method."

    # The Euler step itself.
    x = x0
    y = copy.copy(y0)  # ...so that the user's function cannot mangle y0
    k1 = dydx(x, y, n)
    x1 = x0 + h
    y1 = [0.0] * n
    for i in range(n):
        y1[i] = y0[i] + h * k1[i]
        
    # As a rough error estimate, take a step using a sample
    # point midway and then look at the difference between this
    # and the original Euler step.
    # We will probably underestimate the error as both estimates
    # may have a similar bias.
    err = [0.0] * n
    x2 = x0 + 0.5 * h
    y2 = [0.0] * n
    for i in range(n): y2[i] = 0.5*(y0[i] + y1[i])
    k2 = dydx(x2, y2, n)
    for i in range(n): err[i] = abs(h * (k2[i] - k1[i]))

    # Even though the k2 data should contain better estimates of
    # the y-derivatives, return the original Euler update.
    return (x1, y1, err)


def rkf45_step(x0, h, dydx, n, y0):
    "Steps the set of ODEs by the Runge-Kutta-Fehlberg method."

    # Build up the sample point information.
    x = x0
    y = copy.copy(y0)  # Set up a new list so we don't mangle y0 itself.
    k1 = dydx(x, y, n)
    for i in range(n): k1[i] *= h

    x = x0 + h/4.0
    for i in range(n): y[i] = y0[i] + k1[i]/4.0
    k2 = dydx(x, y, n)
    for i in range(n): k2[i] *= h

    x = x0 + 3.0*h/8.0
    for i in range(n): y[i] = y0[i] + 3.0*k1[i]/32.0 + 9.0*k2[i]/32.0
    k3 = dydx(x, y, n)
    for i in range(n): k3[i] *= h
    
    x = x0 + 12.0*h/13.0
    for i in range(n):
        y[i] = y0[i] + 1932.0*k1[i]/2197.0 - 7200.0*k2[i]/2197.0 \
            + 7296.0*k3[i]/2197.0
    k4 = dydx(x, y, n)
    for i in range(n): k4[i] *= h
    
    x = x0 + h
    for i in range(n):
        y[i] = y0[i] + 439.0*k1[i]/216.0 - 8.0*k2[i] \
            + 3680.0*k3[i]/513.0 - 845.0*k4[i]/4104.0
    k5 = dydx(x, y, n)
    for i in range(n): k5[i] *= h
    
    x = x0 + h/2.0
    for i in range(n):
        y[i] = y0[i] - 8.0*k1[i]/27.0 + 2.0*k2[i] \
            - 3544.0*k3[i]/2565.0 + 1859.0*k4[i]/4104.0 \
            - 11.0*k5[i]/40.0
    k6 = dydx(x, y, n)
    for i in range(n): k6[i] *= h

    # Now, do the integration as a weighting of the sampled data.
    x1 = x0 + h
    y1 = [0.0] * n
    err = [0.0] * n
    for i in range(n):
        y1[i] = y0[i] + 16.0*k1[i]/135.0 + 6656.0*k3[i]/12825.0 \
               + 28561.0*k4[i]/56430.0 - 9.0*k5[i]/50.0 \
               + 2.0*k6[i]/55.0
        err[i] = abs(k1[i]/360.0 - 128.0*k3[i]/4275.0 \
                     - 2197.0*k4[i]/75240.0 + k5[i]/50.0 \
                     + 2.0*k6[i]/55.0)

    return (x1, y1, err)
#----------------------------------------------------------------------

if __name__ == "__main__":
    def dydx_sample_1(x, y, n):
        "Constant derivatives"
        assert n == 3
        return [1.0, -2.0, 3.0]

    def dydx_sample_2(x, y, n):
        "Second-order linear ODE with sine solution"
        assert n == 2
        return [y[1], -y[0]]

    print "Start sample integrations..."
    print "(1) Constant derivatives with Euler method:"
    (x1, y1, err1) = ode_integrate(0.0, 10.0, 1.5,
                                   dydx_sample_1, 3, [0.0]*3, "euler")
    print "x1=", x1
    print "y1=", y1
    print "err1=", err1

    print "(2) Second-order linear ODE with Euler method:"
    (x1, y1, err1) = ode_integrate(0.0, 2.0 * math.pi, 0.01,
                                   dydx_sample_2, 2, [0.0, 1.0], "euler")
    print "x1=", x1
    print "y1=", y1
    print "err1=", err1

    print "(3) Second-order linear ODE with RKF45 method:"
    (x1, y1, err1) = ode_integrate(0.0, 2.0 * math.pi, 0.01,
                                   dydx_sample_2, 2, [0.0, 1.0], "rkf45")
    print "x1=", x1
    print "y1=", y1
    print "err1=", err1
    print "Done."

