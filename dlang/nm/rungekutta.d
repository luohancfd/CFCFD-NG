/**
 * rungekutta.d
 *
 * Ordinary differential equation integration.
 *
 * Author: Peter J.
 * Version: 2014-Jun-15, adapted from the mech2700 class example
 *
 */

module rungekutta;
import std.math;
import std.typecons;

alias Tuple!(double, "t1",  double[], "y1", double[], "err") ODEStepResult;

/**
 * Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
 *
 * Params:
 *     f: a callable function that returns the derivative of y wrt t
 *        The signature of this function is f(t, y) where
 *        t is a float value, y is an array of double values.
 *     t0: is the starting value of the independent variable
 *     y0: an array of starting values for the dependent variables
 *         It is assumed that the y-elements are indexed 0..n-1
 *         where n = y0.length
 *     h: the requested step size
 *
 * Returns:
 * final values of t, y, and error estimates for y values as a tuple.    
 */
ODEStepResult rkf45_step(alias f)(double t0, double[] y0, double h)
    if ( is(typeof(f(0.0, [0.0,0.0])) == double[]) )
{
    // Assuming a system of equations, we need arrays for the intermediate data.
    double[] k1 = y0.dup; double[] k2 = y0.dup; double[] k3 = y0.dup; 
    double[] k4 = y0.dup; double[] k5 = y0.dup; double[] k6 = y0.dup;
    double[] err = y0.dup; double[] y1 = y0.dup;
    // Build up the sample point information as per the text book descriptions.
    k1[] = f(t0, y0.dup);
    k2[] = f(t0 + h/4.0, y0[] + 0.25*h*k1[]);
    k3[] = f(t0 + 3.0*h/8.0, y0[] + 3.0*h*k1[]/32.0 + 9.0*h*k2[]/32.0);
    k4[] = f(t0 + 12.0*h/13.0, y0[] + 1932.0*h*k1[]/2197.0 - 
	     7200.0*h*k2[]/2197.0 + 7296.0*h*k3[]/2197.0);
    k5[] = f(t0 + h, y0[] + 439.0*h*k1[]/216.0 - 8.0*h*k2[] + 
	     3680.0*h*k3[]/513.0 - 845.0*h*k4[]/4104.0);
    k6[] = f(t0 + h/2.0, y0[] - 8.0*h*k1[]/27.0 + 2.0*h*k2[] - 
	     3544.0*h*k3[]/2565.0 + 1859.0*h*k4[]/4104.0 - 11.0*h*k5[]/40.0);
    // Now, do the integration as a weighting of the sampled data.
    double t1 = t0 + h;
    y1[] = y0[] + 16.0*h*k1[]/135.0 + 6656.0*h*k3[]/12825.0 + 
	28561.0*h*k4[]/56430.0 - 9.0*h*k5[]/50.0 + 2.0*h*k6[]/55.0;
    err[] = h*k1[]/360.0 - 128.0*h*k3[]/4275.0 - 2197.0*h*k4[]/75240.0 + 
	h*k5[]/50.0 + 2.0*h*k6[]/55.0;
    foreach(ref e; err) e = fabs(e);
    return ODEStepResult(t1, y1, err);
} // end rkf45_step()

