"""
secant_method.py: Function solver, using the secant method.

Author: PAJ

Version: 17-May-04
"""

def solve(f, x1, x2, tol=1.0e-9):
    """
    Computes x that satisfies f(x) = 0,
    given f and two initial guesses x1 and x2.
    """
    assert callable(f)
    f1 = f(x1)
    f2 = f(x2)
    # Ensure that 2 is closer to the solution than 1
    # and that m3 is the current best guess
    if abs(f2) > abs(f1):
        x3, x1 = x1, x2
        x2 = x3
        f3, f1 = f1, f2
        f2 = f3
    else:
        x3 = x2
        f3 = f2
    # At this point x2 and x3 are equal and
    # should be a better guess than x1.
    # Start the iterative improvement.
    count = 0
    while (abs(f3) > tol) and (count < 100):
        slope = (f2 - f1) / (x2 - x1)
        x3 = x2 - f2 / slope
        f3 = f(x3)
        if abs(f3) < abs(f2):
            # x3 is a better guess than both x2 and x1
            x1, x2 = x2, x3
            f1, f2 = f2, f3
        else:
            # x3 is not as good as x2
            x1 = x3
            f1 = f3
        count = count + 1
    # End of the secant iterations
    if count >= 100:
        print 'secant_method: too many iterations: ', count
        print 'Current guess=', x3, ' Absolute function error=', abs(f3)
    return x3

if __name__ == '__main__':
    print "Begin secant_method test..."
    from math import sin, pi
    def ftest(x): return (sin(x))
    print "one solution at x=", solve(ftest, 6.0, 6.5)
    print "expected x=", (2.0*pi)
    print "Done."
