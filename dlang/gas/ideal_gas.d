/**
 * ideal_gas.d
 * Ideal gas model for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22: initial cut, to explore options.
 */

module ideal_gas;
import gas_model;
import std.math;
import std.stdio;
import std.file;
import std.json;

class Ideal_gas: Gas_model {
public:
    this() {
	// Default model is ideal air, as initialized in the private data.
	_n_species = 1;
	_n_modes = 1;
    }
    this(in char[] file_name) {
	// writefln("Read my parameters from the gas-model file: %s", file_name);
	this();
	if (file_name.length > 0) {
	    auto text = cast(string) read(file_name);
	    auto items = parseJSON(text);
	    _Mmass = items["M"].floating;
	    _gamma = items["gamma"].floating;
	    _Cv = R_universal/_Mmass / (_gamma - 1.0);
	    _Cp = R_universal/_Mmass * _gamma/(_gamma - 1.0);
	    // Reference values for entropy
	    _s1 = items["s1"].floating;
	    _T1 = items["T1"].floating;
	    _p1 = items["p1"].floating;
	    // Molecular transport coefficent constants.
	    _mu_ref = items["mu_ref"].floating;
	    _T_ref = items["T_ref"].floating;
	    _S_mu = items["S_mu"].floating;
	    _k_ref = items["k_ref"].floating;
	    _S_k = items["S_k"].floating;
	}
    }

    override void eval_thermo_state_pT(ref Gas_data Q) {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.rho = Q.p / (Q.T[0] * _Rgas);
	Q.e[0] = Q.T[0] * _Cv;
    }
    override void eval_thermo_state_rhoe(ref Gas_data Q) {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.T[0] = Q.e[0] / _Cv;
	Q.p = Q.rho * _Rgas * Q.T[0];
    }
    override void eval_thermo_state_rhoT(ref Gas_data Q) {
	throw new Exception("not implemented");
    }
    override void eval_thermo_state_rhop(ref Gas_data Q) {
	throw new Exception("not implemented");
    }
    override void eval_thermo_state_ps(ref Gas_data Q, double s) {
	throw new Exception("not implemented");
    }
    override void eval_thermo_state_hs(ref Gas_data Q, double s) {
	throw new Exception("not implemented");
    }
    override void eval_sound_speed(ref Gas_data Q) {
	Q.a = sqrt(_gamma * _Rgas * Q.T[0]);
    }
    override void eval_transport_coefficients(ref Gas_data Q) {
	assert(Q.k.length == 1, "incorrect number of modes");
	Q.mu = _mu_ref * sutherland(Q.T[0], _T_ref, _S_mu);
	Q.k[0] = _k_ref * sutherland(Q.T[0], _T_ref, _S_k);
    }
    override void eval_diffusion_coefficients(ref Gas_data Q) {
	throw new Exception("not implemented");
    }

    override double dedT_const_v(in Gas_data Q) {
	return _Cv;
    }
    override double dhdT_const_p(in Gas_data Q) {
	return _Cp;
    }
    override double gas_constant(in Gas_data Q) {
	return R_universal/_Mmass;
    }
    override double internal_energy(in Gas_data Q) {
	return Q.e[0];
    }
    override double enthalpy(in Gas_data Q) {
	return Q.e[0] + Q.p/Q.rho;
    }
    override double entropy(in Gas_data Q) {
	return _s1 + _Cp * log(Q.T[0]/_T1) - _Rgas * log(Q.p/_p1);
    }

private:
    // Thermodynamic constants
    double _Mmass = 0.02896; // effective molecular mass kg/mole
    double _Rgas = R_universal/0.02896; // J/kg/K
    double _gamma = 1.4;   // ratio of specific heats
    double _Cv = R_universal/0.02896 / 0.4; // J/kg/K
    double _Cp = R_universal/0.02896 * 1.4/0.4; // J/kg/K
    // Reference values for entropy
    double _s1 = 0.0; // J/kg/K
    double _T1 = 298.15; // K
    double _p1 = 101.325e3; // Pa
    // Molecular transport coefficent constants.
    double _mu_ref = 1.716e-5; // Pa.s
    double _T_ref = 273.0; // degrees K
    double _S_mu = 111.0; // degrees K
    double _k_ref = 0.0241; // W/(m.K) 
    double _S_k = 194.0; // degrees K

    /**
     * Sutherland's expression relating the quantity, at temperature T,
     * to the value at the reference temperature.
     *
     * Params:
     *     T: temperature in degrees K
     *     T_ref: reference temperature (degrees K)
     *     S: Sutherland constant (degrees K)
     * Returns:
     *     ratio of property at temperature T to reference value.
     */
    double sutherland(double T, double T_ref, double s) {
	return sqrt(T/T_ref)*(T/T_ref)*(T_ref + s)/(T + s);
    }

} // end class Ideal_gas


unittest {
    import std.stdio;
    auto gm = new Ideal_gas();
    auto gd = new Gas_data(gm);
    assert(approxEqual(gm.R(gd), 287.086), "gas constant");
    assert(gm.n_modes == 1, "number of energy modes");
    assert(gm.n_species == 1, "number of species");
    assert(approxEqual(gd.p, 1.0e5), "pressure");
    assert(approxEqual(gd.T[0], 300.0), "static temperature");
    assert(approxEqual(gd.massf[0], 1.0), "massf[0]");

    gm.eval_thermo_state_pT(gd);
    gm.eval_sound_speed(gd);
    assert(approxEqual(gd.rho, 1.16109), "density");
    assert(approxEqual(gd.e[0], 215314.0), "internal energy");
    assert(approxEqual(gd.a, 347.241), "density");
    gm.eval_transport_coefficients(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), "viscosity");
    assert(approxEqual(gd.k[0], 0.0262449), "conductivity");
}
