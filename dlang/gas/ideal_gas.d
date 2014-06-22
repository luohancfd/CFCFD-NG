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

class Ideal_gas: Gas_model {
public:
    this() {
	// Default model is ideal air
	_n_species = 1;
	_n_modes = 1;
    }
    this(string file_name="gas_model.json") {
	// Read my parameters from the gas-model file.
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
    double _Mmass = 28.96; // effective molecular mass
    double _Rgas = R_universal/28.96; // J/kg/K
    double _gamma = 1.4;   // ratio of specific heats
    double _Cv = R_universal/28.96 / 0.4; // J/kg/K
    double _Cp = R_universal/28.96 * 1.4/0.4; // J/kg/K
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
	return pow(T/T_ref, 1.5) * (T_ref + s)/(T + s);
    }

} // end class Ideal_gas
