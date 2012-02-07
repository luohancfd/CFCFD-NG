"""
sode.py: Integrate a set of stiff ODEs.

Author: P.Jacobs

Version: 1.0, 21-Feb-2005

Transcript::

    Start sample integrations...
    (1) Constant derivatives with semi-implicit RK1 method:
    t1= 10.0
    y1= [ 10. -20.  30.]
    err1= [ 0.  0.  0.]
    (2) Second-order linear ODE with semi-implicit RK1 method:
    t1= 6.28318530718
    y1= [ -5.23352415e-05   9.99999999e-01]
    err1= [ 0.01999997  0.01998885]
    (3) Second-order linear ODE with semi-implicit RK2 method:
    t1= 6.28318530718
    y1= [ -2.29779805e-10   1.00000000e-00]
    err1= [ 0.01999992  0.01998906]
    Done.
"""

from math import sqrt, pi

try:
    from numpy import array, identity
    from numpy import float as Float
    from numpy.linalg import solve as solve_linear_equations
except:
    try:
        from Numeric import array, identity, Float
        from LinearAlgebra import solve_linear_equations
    except:
        print "Could import neither numpy nor Numeric."


def ode_integrate(t0, tlast, dt, F, n, y0, dFdt=None, dFdy=None, method="rk1"):
    """
    Steps the set of ODEs until x reaches xlast.

    This function coordinates the work of integrating a system
    of first-order differential equations of the form:

        | y' = F(t, y)
        | y(t=t0) = y0

    The actual work is done by one of a set of more specialised
    stepping functions that appear below.

    :param t0: the starting value of the independent variable
    :param tlast: the desired finishing value for x
    :param dt: the requested step size
    :param F: a callable function that returns the derivative of y wrt t
              The signature of this function is F(t, y, n) where
              t is a float value, y is a list (or array) or float values
              and n is an integer specifying the number of equations.
    :param n: the number of dependent variables (in y)
    :param y0: an array of starting values for the dependent variables
               It is assumed that the y-elements are indexed 0..n-1
    :param method: a string specifying which stepping method to use.
                   "rk1", "rk2"

    :returns: final values of t, y, and error estimates for y values as a tuple.    
    """
    assert callable(F)
    assert n <= len(y0)
    assert t0 <= tlast
    t = t0
    y = y0.copy()  # Set up a new list so we don't mangle y0 itself.
    err_sum = array([0.0] * n)
    while t < tlast:
        # Decide if we need to take a full step
        if t + dt <= tlast:
            h = dt
        else:
            h = tlast - t
        # The real work is delegated.    
        if method == "rk1":
            (t, y, err) = rk1_step(t, h, F, n, y, dFdt, dFdy)
        elif method == "rk2":
            (t, y, err) = rk2_step(t, h, F, n, y, dFdt, dFdy)
        else:
            raise ValueError, "Requested stepping method is not available."
        err_sum += err
    # Return the point reached, the set of y-values and an error estimate.
    return (t, y, err_sum)


def rk1_step(t0, h, F, n, y0, dFdt, dFdy):
    """
    Steps the set of ODEs by the semi-implicit one-stage Runge-Kutta method.
    """
    t = t0
    y = y0.copy()  # ...so that the user's function cannot mangle y0
    k0 = F(t, y, n)
    rhs = k0 + dFdt(t, y, n) * h / 2.0
    lhs = identity(n, Float) - dFdy(t, y, n) * h / 2.0
    k1 = solve_linear_equations(lhs, rhs)
    t1 = t0 + h
    y1 = y0 + h * k1
    # As a poor error estimate, look at the difference between this rk1 step
    # and an explicit-Euler step.
    err = abs(h * (k0 - k1))
    return (t1, y1, err)


def rk2_step(t0, h, F, n, y0, dFdt, dFdy):
    """
    Steps the set of ODEs by the semi-implicit two-stage Runge-Kutta method.
    """
    t = t0
    y = y0.copy()  # ...so that the user's function cannot mangle y0
    k0 = F(t, y, n)
    k1 = k0.copy(); k1_old = k0.copy()
    k2 = k0.copy(); k2_old = k0.copy()
    # constants used in computing sample points
    c1t = 0.5*(1 - 1.0/sqrt(3.0))
    c2t = 0.5*(1 + 1.0/sqrt(3.0))
    c1y = 0.5*(0.5 - 1.0/sqrt(3.0))
    c2y = 0.5*(0.5 + 1.0/sqrt(3.0))
    # Apply fixed-point iteration to get k1, k2.
    converged = False
    count = 0
    TOL = 1.0e-7
    while not converged and count < 100:
        k1 = F(t0 + c1t*h, y0 + h*(0.25*k1 + c1y*k2), n)
        k2 = F(t0 + c2t*h, y0 + h*(c2y*k1 + 0.25*k2), n)
        d1 = max(abs((k1 - k1_old)/(k1 + 1.0)))
        d2 = max(abs((k2 - k2_old)/(k2 + 1.0)))
        converged = d1 < TOL and d2 < TOL
        count += 1
        k1_old = k1.copy()
        k2_old = k2.copy()
    if count >= 100:
        print "rk2_step(): Warning: iteration for k1, k2 did not converge."
    t1 = t0 + h
    y1 = y0 + h*0.5*(k1 + k2)
    # As a poor error estimate, look at the difference between this rk2 step
    # and an explicit-Euler step.
    err = abs(h * (k0 - 0.5*(k1 + k2)))
    return (t1, y1, err)

#--------------------------------------------------------------------------------

if __name__ == "__main__":
    def F1(t, y, n):
        "Constant derivatives"
        assert n == 3, "n not 3"
        return array([1.0, -2.0, 3.0])

    def dF1dt(t, y, n):
        "List of derivatives."
        return array([0.0, 0.0, 0.0])

    def dF1dy(t, y, n):
        "Jacobian"
        return array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    
    def F2(t, y, n):
        "Second-order linear ODE with sine solution"
        assert n == 2, "n not 2"
        return array([y[1], -y[0]])

    def dF2dt(t, y, n):
        "List of derivatives."
        return array([0.0, 0.0])

    def dF2dy(t, y, n):
        "Jacobian"
        return array([[0.0, 1.0], [-1.0, 0.0]])

    print "Start sample integrations..."
    print "(1) Constant derivatives with semi-implicit RK1 method:"
    (t1, y1, err1) = ode_integrate(0.0, 10.0, 1.5, F1, 3, array([0.0]*3),
                                   dF1dt, dF1dy, "rk1")
    print "t1=", t1
    print "y1=", y1
    print "err1=", err1

    print "(2) Second-order linear ODE with semi-implicit RK1 method:"
    (t1, y1, err1) = ode_integrate(0.0, 2.0*pi, 0.01, F2, 2,
                                   array([0.0, 1.0]), dF2dt, dF2dy, "rk1")
    print "t1=", t1
    print "y1=", y1
    print "err1=", err1

    print "(3) Second-order linear ODE with semi-implicit RK2 method:"
    (t1, y1, err1) = ode_integrate(0.0, 2.0*pi, 0.01, F2, 2,
                                   array([0.0, 1.0]), dF2dt, dF2dy, "rk2")
    print "t1=", t1
    print "y1=", y1
    print "err1=", err1

    print "Done."
