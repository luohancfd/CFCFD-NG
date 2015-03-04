/**
 * reaction_mechanism.d
 *
 * Author: Rowan G.
 */

module kinetics.reaction_mechanism;

import std.stdio;

import gas;
import luad.all;
import util.lua_service;
import util.msg_service;
import kinetics.rate_constant;
import kinetics.reaction;

class ReactionMechanism
{
public:
    @property ulong n_reactions() const { return _reactions.length; }
    this(in Reaction[] reactions, size_t n_species)
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
    ReactionMechanism dup() const
    {
	return new ReactionMechanism(_reactions, _q.length);
    }
    
    final void eval_rate_constants(in GasState Q)
    {
	foreach ( ref r; _reactions ) r.eval_rate_constants(Q);
    }

    final void eval_rates(in double[] conc, double[] rates)
    {
	eval_split_rates(conc, _q, _L);
	foreach ( isp; 0..conc.length ) rates[isp] = _q[isp] - _L[isp];
    }
    final void eval_split_rates(in double[] conc, double[] q, double[] L)
    {
	foreach ( ref r; _reactions ) r.eval_rates(conc);
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
    final double k_f(int ir)
    {
	return _reactions[ir].k_f;
    }
    final double k_b(int ir)
    {
	return _reactions[ir].k_b;
    }
    
private:
    Reaction[] _reactions;
    // Working array space
    static bool _work_arrays_initialised = false;
    static double[] _q;
    static double[] _L;
}

/++
 + Creates a ReactionMechanism object from information contained in a LuaTable.
 +
 + The expected format of the table is a of Reaction objects, eg.
 + tab = { [1]={... reaction 1 info ...},
 +         [2]={... reaction 2 info ...},
 +              .......................
 +         [n]={... reaction n info ...}
 +       }
 +/ 
ReactionMechanism createReactionMechanism(LuaTable t, in GasModel gmodel)
{
    auto n_species = gmodel.n_species;
    auto n_reactions = t.length;
    Reaction[] reactions;
    foreach ( i; 1..n_reactions+1 ) {
	reactions ~= createReaction(t.get!LuaTable(i), n_species);
    }
    return new ReactionMechanism(reactions, n_species);
}

ReactionMechanism createReactionMechanism(string fname, in GasModel gmodel)
{
    auto lua = initLuaState(fname);
    return createReactionMechanism(lua.get!LuaTable("reaction"), gmodel);
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
    reacMech.eval_rate_constants(gd);
    reacMech.eval_rates(conc, rates);
    assert(approxEqual([-643.9303, -643.9303, 1287.8606], rates), failedUnitTest());
}
