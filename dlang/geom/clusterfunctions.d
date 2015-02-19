/** clusterfunctions.d
 * Clustering functions ported from the C++ code.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2015-02-20 first code
 */

module clusterfunctions;

import std.math;

double[] distribute_parameter_values(int n,
				     double t0, double t1,
				     double delegate(double) f)
// Returns an array of parameter values, distributed from t0 though t1.
// f determines the form of distribution.
{
    double[] tv;
    double dt = 1.0 / (n-1);
    foreach (i; 0 .. n) {
	double t = f(dt * i); // first mapping is done by the cluster function
	tv ~= (1.0 - t) * t0 + t * t1; // map to specified range and save
    }
    return tv;
}

double delegate(double) linearFunction(double t0, double t1)
// Returns a function that linearly interpolates end values.
{
    auto f = delegate double(double t) {
	return (1.0-t) * t0 + t * t1;
    };
    return f;
}

double delegate(double) robertsFunction(bool end0, bool end1, double beta)
// Returns a function that can be used for transforming the parameter t
// that is used to traverse each dimension for the parametric geometry objects.
{
    // Decide on stretching parameters for original Robert's transform.
    double alpha = 0.5;
    bool reverse = false;
    bool cluster = true;
    if (!end0 && !end1) cluster = false;
    if (beta <= 1.0) cluster = false;
    if (end0 && end1) alpha = 0.5;
    if (end0 && !end1) {
	reverse = true;
	alpha   = 0.0;
    }
    if (!end0 && end1) {
	reverse = 0;
	alpha   = 0.0;
    }
    // Build the actual function that will do the work later.
    auto f = delegate double(double t) {
	double tbar;
	if (reverse) t = 1.0 - t;
	if (cluster) {
	    tbar = roberts_original(t, alpha, beta);
	} else { 
	    tbar = t;
	}
	if (reverse) tbar = 1.0 - tbar;
	return tbar;
    };
    return f;
} // end robertsClusterFunction()

double roberts_original(double eta, double alpha, double beta)
// eta   : unstretched coordinate, 0 <= eta <= 1.0
// alpha : location of stretching
//         alpha = 0.5 : clustering of nodes at both extremes of eta
//         alpha = 0.0 : nodes will be clustered near eta=1.0
// beta  : stretching factor (more stretching as beta --> 1.0)
// Returns stretched coordinate, 0 <= roberts <= 1.0
{
    double lambda = (beta + 1.0) / (beta - 1.0);
    lambda = pow ( lambda, ((eta - alpha)/(1.0 - alpha)) );
    double etabar = (beta + 2.0 * alpha) * lambda - beta + 2.0 * alpha;
    etabar = etabar / ((2.0 * alpha + 1.0) * (1.0 + lambda));
    return etabar;
} // end roberts_original()


unittest {
    auto cf = robertsFunction(false, true, 1.1);
    assert(approxEqual(cf(0.1), 0.166167), "robertsFunction");
    assert(approxEqual(cf(0.9), 0.96657), "robertsFunction");
    auto cf2 = linearFunction(1.0, 0.0);
    assert(approxEqual(cf2(0.1), 0.9), "linearFunction");
    assert(approxEqual(cf2(0.9), 0.1), "linearFunction");
}
