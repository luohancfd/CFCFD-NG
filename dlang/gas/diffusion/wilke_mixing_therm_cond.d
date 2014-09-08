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
import therm_cond;

class WilkeMixingThermCond : ThermalConductivity {
public:
    this(in ThermalConductivity[] tcms, in double[] MW) {
	assert(tcms.length == MW.length);
	foreach (tcm; tcms) {
	    _tcms ~= tcm.dup;
	}
	_MW = MW.dup;
	_x.length = _MW.length;
	_k.length = _MW.length;
	_psi.length = _MW.length;
	foreach (ref p; _psi) {
	    p.length = _MW.length;
	}
    }
    this(in WilkeMixingThermCond src) {
	foreach (tcm; src._tcms) {
	    _tcms ~= tcm.dup;
	}
	_MW = src._MW.dup;
	_x = src._x.dup;
	_k = src._k.dup;
	_psi.length = src._psi.length;
	foreach (i; 0 .. src._psi.length) {
	    _psi[i] = src._psi[i].dup;
	}
    }
    override WilkeMixingThermCond dup() const {
	return new WilkeMixingThermCond(this);
    }

    override void update_thermal_conductivity(ref GasState Q) {
	// 1. Evaluate the mole fractions
	massf2molef(Q.massf, _MW, _x);
	// 2. Calculate the component viscosities
	for ( auto i = 0; i < Q.massf.length; ++i ) {
	    _tcms[i].update_thermal_conductivity(Q);
	    _k[i] = Q.k[0];
	}
	// 3. Calculate interaction potentials
	for ( auto i = 0; i < Q.massf.length; ++i ) {
	    for ( auto j = 0; j < Q.massf.length; ++j ) {
		double numer = pow((1.0 + sqrt(_k[i]/_k[j])*pow(_MW[j]/_MW[i], 0.25)), 2.0);
		double denom = (4.0/sqrt(2.0))*sqrt(1.0 + (_MW[i]/_MW[j]));
		_psi[i][j] = numer/denom;
	    }
	}
	// 4. Apply mixing formula
	double sum;
	double k = 0.0;
	for ( auto i = 0; i < Q.massf.length; ++i ) {
	    if ( _x[i] < SMALL_MOLE_FRACTION ) continue;
	    sum = 0.0;
	    for ( auto j = 0; j < Q.massf.length; ++j ) {
		if ( _x[j] < SMALL_MOLE_FRACTION ) continue;
		sum += _x[j]*_psi[i][j];
	    }
	    k += _k[i]/(1.0 + (1.0/_x[i])*sum);
	}
	Q.k[0] = k;
    }

private:
    ThermalConductivity[] _tcms; // component viscosity models
    double[] _MW; // component molecular weights
    // Working array space
    double[] _x;
    double[] _k;
    double[][] _psi;
}

unittest {
    import std.stdio;
    import sutherland_therm_cond;
    // Placeholder test. Redo with CEA curves.
    double T = 300.0;
    auto tcm_N2 = new SutherlandThermCond(273.0, 0.0242, 150.0);
    auto tcm_O2 = new SutherlandThermCond(273.0, 0.0244, 240.0);
    auto tcm = new WilkeMixingThermCond([tcm_N2, tcm_O2], [28.0e-3, 32.0e-3]);

    auto gd = GasState(2, 1);
    gd.T[0] = T;
    gd.massf[0] = 0.8;
    gd.massf[1] = 0.2;
    tcm.update_thermal_conductivity(gd);
    assert(approxEqual(0.0159736, gd.k[0]));
}
