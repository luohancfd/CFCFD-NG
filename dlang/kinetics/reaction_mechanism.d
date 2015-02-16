/**
 * reaction_mechanism.d
 *
 * Author: Rowan G.
 */

class ReactionMechanism
{
public:
    @property int n_reactions() const { return _reactions.length; }
    
    final void eval_rates(in GasState Q, in double[] conc, bool evalRateConstants, out double[] rate) const
    {
	eval_split_rates(conc, evalRateConstants, _q, _L);
	for ( auto isp; 0..conc.length ) {
	    rate[isp] = _q[isp] - _L[isp];
	}
    }
    final void eval_split_rates(in GasState Q, in double[] conc, bool evalRateConstants,
				out double[] q, out double[] L) const
    {
	if ( evalRateConstants ) {
	    for ( r; _reactions )
		r.eval_rate_constants(Q);
	}
	
	for ( r; _reactions )
	    r.eval_rates(Q);

	for ( isp; 0.._n_species ) {
	    q[isp] = 0.0;
	    L[isp] = 0.0;
	}
	
	for ( ir, r; _reactions ) {
	    for ( isp; _participators[ir] ) {
		q[isp] += r.production(isp);
		L[isp] += r.loss(isp);
	    }
	}
    }

private:
    Reaction[] _reactions;
    // Working array space
    static double[] _q;
    static double[] _L;
}

