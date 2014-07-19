/**
 * idealgas.d
 * Ideal gas model for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22: initial cut, to explore options.
 */

module idealgas;
import gasmodel;
import std.math;
import std.stdio;
import std.file;
import std.json;
import std.conv;

class IdealGas: GasModel {
public:
    this() {
	// Default model is mostly initialized in the private data below.
	_n_species = 1;
	_n_modes = 1;
	_species_names ~= "ideal air";
    }
    this(in char[] file_name) {
	// writefln("Read my parameters from the gas-model file: %s", file_name);
	this();
	if (file_name.length > 0) {
	    auto text = cast(string) read(file_name);
	    auto items = parseJSON(text);
	    _species_names ~= items["name"].str;
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

    override string toString()
    {
	char[] repr;
	repr ~= "IdealGas(";
	repr ~= "name=\"" ~ _species_names[0] ~"\"";
	repr ~= ", Mmass=" ~ to!string(_Mmass);
	repr ~= ", gamma=" ~ to!string(_gamma);
	repr ~= ", s1=" ~ to!string(_s1);
	repr ~= ", T1=" ~ to!string(_T1);
	repr ~= ", p1=" ~ to!string(_p1);
	repr ~= ", mu_ref=" ~ to!string(_mu_ref);
	repr ~= ", T_ref=" ~ to!string(_T_ref);
	repr ~= ", S_mu=" ~ to!string(_S_mu);
	repr ~= ", k_ref=" ~ to!string(_k_ref);
	repr ~= ", S_k=" ~ to!string(_S_k);
	repr ~= ")";
	return to!string(repr);
    }

    override const void update_thermo_from_pT(ref GasState Q) {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.rho = Q.p / (Q.T[0] * _Rgas);
	Q.e[0] = Q.T[0] * _Cv;
    }
    override const void update_thermo_from_rhoe(ref GasState Q) {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.T[0] = Q.e[0] / _Cv;
	Q.p = Q.rho * _Rgas * Q.T[0];
    }
    override const void update_thermo_from_rhoT(ref GasState Q) {
	throw new Exception("not implemented");
    }
    override const void update_thermo_from_rhop(ref GasState Q) {
	throw new Exception("not implemented");
    }
    override const void update_thermo_from_ps(ref GasState Q, double s) {
	throw new Exception("not implemented");
    }
    override const void update_thermo_from_hs(ref GasState Q, double s) {
	throw new Exception("not implemented");
    }
    override const void update_sound_speed(ref GasState Q) {
	Q.a = sqrt(_gamma * _Rgas * Q.T[0]);
    }
    override const void update_trans_coeffs(ref GasState Q) {
	assert(Q.k.length == 1, "incorrect number of modes");
	Q.mu = _mu_ref * sutherland(Q.T[0], _T_ref, _S_mu);
	Q.k[0] = _k_ref * sutherland(Q.T[0], _T_ref, _S_k);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
	throw new Exception("not implemented");
    }
    */
    override const double dedT_const_v(in GasState Q) {
	return _Cv;
    }
    override const double dhdT_const_p(in GasState Q) {
	return _Cp;
    }
    override const double gas_constant(in GasState Q) {
	return R_universal/_Mmass;
    }
    override const double internal_energy(in GasState Q) {
	return Q.e[0];
    }
    override const double enthalpy(in GasState Q) {
	return Q.e[0] + Q.p/Q.rho;
    }
    override const double entropy(in GasState Q) {
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
    const double sutherland(double T, double T_ref, double s) {
	return sqrt(T/T_ref)*(T/T_ref)*(T_ref + s)/(T + s);
    }

} // end class Ideal_gas


unittest {
    import std.stdio;
    auto gm = new IdealGas();
    assert(gm.species_name(0) == "ideal air", "species name list");
    auto gd = new GasState(gm, 100.0e3, 300.0);
    assert(approxEqual(gm.R(gd), 287.086), "gas constant");
    assert(gm.n_modes == 1, "number of energy modes");
    assert(gm.n_species == 1, "number of species");
    assert(approxEqual(gd.p, 1.0e5), "pressure");
    assert(approxEqual(gd.T[0], 300.0), "static temperature");
    assert(approxEqual(gd.massf[0], 1.0), "massf[0]");

    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    assert(approxEqual(gd.rho, 1.16109), "density");
    assert(approxEqual(gd.e[0], 215314.0), "internal energy");
    assert(approxEqual(gd.a, 347.241), "density");
    gm.update_trans_coeffs(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), "viscosity");
    assert(approxEqual(gd.k[0], 0.0262449), "conductivity");
}
