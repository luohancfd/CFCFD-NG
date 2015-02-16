/**
 * reaction.d
 * This module holds the classes related to computing
 * rates of species change due to a single reaction.
 *
 * Author: Rowan J. Gollan
 */

module kinetics.reaction;

import std.math;
import std.typecons;

import gas;
/++
 Reaction is an interface specifying the public services
 provided by an object of Reaction type.
+/
interface Reaction
{
public:
    @property double k_f() const { return _k_f; }
    @property double k_b() const { return _k_b; }
    
    final void eval_rate_constants(in GasState Q)
    {
	_k_f = eval_forward_rate_constant(Q);
	_k_b = eval_backward_rate_constant(Q);
    }

    final void eval_rates(in double[] conc)
    {
	_w_f = eval_forward_rate(conc);
	_w_b = eval_backward_rate(conc);
    }

    double production(int isp) const;
    double loss(int isp) const;
	    
private:
    RateConstant _forward, _backward;
    double _k_f, _k_b; // Storage of computed rate constants
    double _w_f, _w_b; // Storage of computed rates of change
    double eval_forward_rate_constant(in GasState Q) const
    {
	return _forward.eval(Q);
    }
    double eval_backward_rate_constant(in GasState Q) const
    {
	return _backward.eval(Q);
    }
    double eval_forward_rate(in double[] conc) const;
    double eval_backward_rate(in double[] conc) const;
}

/++
 An ElementaryReaction is a reaction whose rate law can
 be written directly from its molecularity.
 There are three forms:
   - unimolecular
   - bimolecular
   - termolecular
+/
class ElementaryReaction : Reaction
{
public:
    override double production(isp) const
    {
	if ( _nu[isp] >  0 )
	    return _nu[isp]*_w_f;
	else if ( _nu[isp] < 0 )
	    return -_nu[isp]*_w_b;
	else
	    return 0.0;
    }
    
    override double loss(isp) const
    {
	if ( _nu[isp] > 0 ) 
	    return _nu[isp]*_w_b;
	else if ( _nu[isp] < 0 )
	    return -_nu[isp]*_w_f;
	else 
	    return 0.0
    }

private:
    Tuple!(int, int)[] _reactants;
    Tuple!(int, int)[] _products;
    int[] _nu;
    override double eval_forward_rate(in double[] conc) const
    {
	double val = _k_f;
	for ( r; _reactants ) {
	    int isp = r[0];
	    int coeff = r[1];
	    val *= pow(conc[isp], coeff);
	}
	return val;
    }

    override double eval_backward_rate(in double[] conc) const
    {
	double val = _k_b;
	for ( p; _products ) {
	    int isp = p[0];
	    int coeff = p[1];
	    val *= pow(conc[isp], coeff);
	}
	return val;
    }
}
