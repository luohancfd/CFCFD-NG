/**
 * wilke_mixing.d
 * Implements the Gordon and McBride variant
 * of Wilke's mixing rule to compute the
 * viscosity and thermal conductivity of
 * a mixture of gases.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2014-09-08 -- initial cut
 */

import std.math;
import gasmodel;
import viscosity;

class WilkeMixingViscosity : Viscosity {
public:
    this(in Viscosity[] vms, in double[] MW) {
	assert(vms.length == MW.length);
	foreach (v; vms) {
	    _vms ~= v.dup;
	}
	_MW = MW.dup;
	_x.length = _MW.length;
	_mu.length = _MW.length;
	_phi.length = _MW.length;
	foreach (ref p; _phi) {
	    p.length = _MW.length;
	}
    }
    this(in WilkeMixingViscosity src) {
	foreach (v; src._vms) {
	    _vms ~= v.dup;
	}
	_MW = src._MW.dup;
	_x = src._x.dup;
	_mu = src._mu.dup;
	_phi.length = src._phi.length;
	foreach (i; 0 .. src._phi.length) {
	    _phi[i] = src._phi[i].dup;
	}
    }
    override WilkeMixingViscosity dup() const {
	return new WilkeMixingViscosity(this);
    }

    override void update_viscosity(ref GasState Q) {
	// 1. Evaluate the mole fractions
	massf2molef(Q.massf, _MW, _x);
	// 2. Calculate the component viscosities
	for ( auto i = 0; i < Q.massf.length; ++i ) {
	    _vms[i].update_viscosity(Q);
	    _mu[i] = Q.mu;
	}
	// 3. Calculate interaction potentials
	for ( auto i = 0; i < Q.massf.length; ++i ) {
	    for ( auto j = 0; j < Q.massf.length; ++j ) {
		double numer = pow((1.0 + sqrt(_mu[i]/_mu[j])*pow(_MW[j]/_MW[i], 0.25)), 2.0);
		double denom = (4.0/sqrt(2.0))*sqrt(1.0 + (_MW[i]/_MW[j]));
		_phi[i][j] = numer/denom;
	    }
	}
	// 4. Apply mixing formula
	double sum;
	double mu = 0.0;
	for ( auto i = 0; i < Q.massf.length; ++i ) {
	    if ( _x[i] < SMALL_MOLE_FRACTION ) continue;
	    sum = 0.0;
	    for ( auto j = 0; j < Q.massf.length; ++j ) {
		if ( _x[j] < SMALL_MOLE_FRACTION ) continue;
		sum += _x[j]*_phi[i][j];
	    }
	    mu += _mu[i]/(1.0 + (1.0/_x[i])*sum);
	}
	Q.mu = mu;
    }

private:
    Viscosity[] _vms; // component viscosity models
    double[] _MW; // component molecular weights
    // Working array space
    double[] _x;
    double[] _mu;
    double[][] _phi;
}

unittest {
    import sutherland_visc;
    // Placeholder test. Redo with CEA curves.
    double T = 300.0;
    auto vm_N2 = new SutherlandViscosity(273.0, 1.663e-5, 107.0);
    auto vm_O2 = new SutherlandViscosity(273.0, 1.919e-5, 139.0);
    auto vm = new WilkeMixingViscosity([vm_N2, vm_O2], [28.0e-3, 32.0e-3]);

    auto gd = GasState(2, 1);
    gd.T[0] = T;
    gd.massf[0] = 0.8;
    gd.massf[1] = 0.2;
    vm.update_viscosity(gd);
    assert(approxEqual(1.12102e-05, gd.mu));
}
