## \file ode.py
"""\brief Integrate a set of first-order ODEs.

This module provides a small number of ODE integration schemes
together with a coordination function ode_integrate().
The schemes will integrate a system of (nonstiff) ODEs
over a given range of the independent variable and with
a specified step.
The user is expected to supply a function that accepts the
current value of the independent variable and system state
and which provides the derivatives of the system as an array.
For further details, see the doc string of ode_integrate() and study
the sample functions near the end of this file.

Running this module as a Python script gives me the following transcript:

Start sample integrations...
(1) Constant derivatives with Euler method:
t1= 10.0
y1= [ 10. -20.  30.]
err1= [ 0.  0.  0.]
(2) Second-order linear ODE with Euler method:
t1= 6.28318530718
y1= [ -2.16009598e-04   1.03190177e+00]
err1= [ 0.02031786  0.02030512]
(3) Second-order linear ODE with RKF45 method:
t1= 6.28318530718
y1= [  8.42476994e-14   1.00000000e+00]
err1= [  5.12421830e-11   5.12821769e-11]
Done.

\author P. Jacobs
        School of Engineering, UQ
\version 05-Oct-03 implementation with lists for storage
\version 21-Feb-05 use Numeric arrays for storage and manipulation
                   of the state data.
"""

from math import pi
try:
    from numpy import array
except:
    try:
        from Numeric import array
    except:
        print "Could import neither numpy nor Numeric."


def ode_integrate(t0, tlast, dt, F, n, y0, method="euler"):
    """Steps the set of ODEs until t reaches tlast.

    This function coordinates the work of integrating a system
    of first-order differential equations of the form
        y' = f(t, y)
        y(t=t0) = y0
    The actual work is done by one of a set of more specialised
    stepping functions that appear below.

    \param t0    is the starting value of the independent variable
    \param tlast is the desired finishing value for x
    \param dt    is the requested step size
    \param F     is a callable function that returns the derivative of y wrt t
                 The signature of this function is F(t, y, n) where
                 t is a float value, y is a list (or array) or float values
                 and n is an integer specifying the number of equations.
    \param n     is the number of dependent variables (in y)
    \param y0    is an array of starting values for the dependent variables
                 It is assumed that the y-elements are indexed 0..n-1
    \param method is a string specifying which stepping method to use.

    \returns final values of t, y, and error estimates for y values
             are returned as a tuple.    
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
        if method == "euler":
            (t, y, err) = euler_step(t, h, F, n, y)
        elif method == "rkf45":
            (t, y, err) = rkf45_step(t, h, F, n, y)
        else:
            raise ValueError, "Requested stepping method is not available."
        err_sum += err
    # Return the point reached, the set of y-values and an error estimate.
    return (t, y, err_sum)


def euler_step(t0, h, F, n, y0):
    "Steps the set of ODEs by the Euler method."
    # The Euler step itself.
    t = t0
    y = y0.copy()  # ...so that the user's function cannot mangle y0
    k1 = F(t, y, n)
    t1 = t0 + h
    y1 = y0 + h * k1
    # As a rough error estimate, take a step using a sample
    # point midway and then look at the difference between this
    # and the original Euler step.
    # We will probably underestimate the error as both estimates
    # may have a similar bias.
    t2 = t0 + 0.5 * h
    y2 = 0.5*(y0 + y1)
    k2 = F(t2, y2, n)
    err = abs(h * (k2 - k1))
    # Even though the k2 data should contain better estimates of
    # the y-derivatives, return the original Euler update.
    return (t1, y1, err)


def rkf45_step(t0, h, F, n, y0):
    "Steps the set of ODEs by the Runge-Kutta-Fehlberg method."
    # Build up the sample point information.
    k1 = F(t0, y0.copy(), n)
    k2 = F(t0 + h/4.0, y0 + 0.25*h*k1, n)
    k3 = F(t0 + 3.0*h/8.0, y0 + 3.0*h*k1/32.0 + 9.0*h*k2/32.0, n)
    k4 = F(t0 + 12.0*h/13.0, y0 + 1932.0*h*k1/2197.0 - 7200.0*h*k2/2197.0 + 7296.0*h*k3/2197.0, n)
    k5 = F(t0 + h, y0 + 439.0*h*k1/216.0 - 8.0*h*k2 + 3680.0*h*k3/513.0 - 845.0*h*k4/4104.0, n)
    k6 = F(t0 + h/2.0, y0 - 8.0*h*k1/27.0 + 2.0*h*k2 - 3544.0*h*k3/2565.0 + 1859.0*h*k4/4104.0 - 11.0*h*k5/40.0, n)
    # Now, do the integration as a weighting of the sampled data.
    t1 = t0 + h
    y1 = y0 + 16.0*h*k1/135.0 + 6656.0*h*k3/12825.0 + 28561.0*h*k4/56430.0 - 9.0*h*k5/50.0 + 2.0*h*k6/55.0
    err = abs(h*k1/360.0 - 128.0*h*k3/4275.0 - 2197.0*h*k4/75240.0 + h*k5/50.0 + 2.0*h*k6/55.0)
    return (t1, y1, err)

#----------------------------------------------------------------------

if __name__ == "__main__":
    def F_sample_1(t, y, n):
        "Constant derivatives"
        assert n == 3
        return array([1.0, -2.0, 3.0])

    def F_sample_2(t, y, n):
        "Second-order linear ODE with sine solution"
        assert n == 2
        return array([y[1], -y[0]])

    print "Start sample integrations..."
    print "(1) Constant derivatives with Euler method:"
    (t1, y1, err1) = ode_integrate(0.0, 10.0, 1.5, F_sample_1, 3, array([0.0]*3), "euler")
    print "t1=", t1
    print "y1=", y1
    print "err1=", err1

    print "(2) Second-order linear ODE with Euler method:"
    (t1, y1, err1) = ode_integrate(0.0, 2.0*pi, 0.01, F_sample_2, 2, array([0.0, 1.0]), "euler")
    print "t1=", t1
    print "y1=", y1
    print "err1=", err1

    print "(3) Second-order linear ODE with RKF45 method:"
    (t1, y1, err1) = ode_integrate(0.0, 2.0*pi, 0.01, F_sample_2, 2, array([0.0, 1.0]), "rkf45")
    print "t1=", t1
    print "y1=", y1
    print "err1=", err1
    print "Done."

