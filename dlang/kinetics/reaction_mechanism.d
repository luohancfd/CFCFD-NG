/**
 * reaction_mechanism.d
 *
 * Author: Rowan G.
 */

module kinetics.reaction_mechanism;

import gas;
import util.msg_service;
import kinetics.rate_constant;
import kinetics.reaction;

class ReactionMechanism
{
public:
    @property ulong n_reactions() const { return _reactions.length; }
    this(in Reaction[] reactions, int n_species)
    {
	foreach ( ref r; reactions) {
	    _reactions ~= r.dup();
	}
	// Allocate work arrays
	if (!_work_arrays_initialised) {
	    _q.length = n_species;
	    _L.length = n_species;
	    _work_arrays_initialised = true;
	}
    }
    
    final void eval_rates(in GasState Q, in double[] conc, bool evalRateConstants, out double[] rate)
    {
	eval_split_rates(Q, conc, evalRateConstants, _q, _L);
	foreach ( isp; 0..conc.length ) {
	    rate[isp] = _q[isp] - _L[isp];
	}
    }
    final void eval_split_rates(in GasState Q, in double[] conc, bool evalRateConstants,
				out double[] q, out double[] L)
    {
	if ( evalRateConstants ) {
	    foreach ( ref r; _reactions )
		r.eval_rate_constants(Q);
	}
	foreach ( ref r; _reactions )
	    r.eval_rates(conc);
	foreach ( isp; 0..conc.length ) {
	    q[isp] = 0.0;
	    L[isp] = 0.0;
	}
	foreach ( ir, ref r; _reactions ) {
	    foreach ( isp; r.participants ) {
		q[isp] += r.production(isp);
		L[isp] += r.loss(isp);
	    }
	}
    }

private:
    Reaction[] _reactions;
    // Working array space
    static bool _work_arrays_initialised = false;
    static double[] _q;
    static double[] _L;
}

unittest
{
    import std.math;
    // Test the rate of concentration change at the initial
    // condition for the H2 + I2 reaction system.
    double[] conc = [4.54, 4.54, 0.0];
    auto rc = new ArrheniusRateConstant(1.94e14, 0.0, 20620.0);
    auto gd = new GasState(3, 1);
    gd.T[0] = 700.0;
    auto reaction = new ElementaryReaction(rc, rc, [0, 1], [1, 1],
					   [2], [2], 3);
    auto reacMech = new ReactionMechanism([reaction], 3);
    double[] rates;
    rates.length = 3;
    reacMech.eval_rates(gd, conc, true, rates);
    assert(approxEqual([-643.9303, -643.9303, 1287.8606], rates), failedUnitTest());
}
