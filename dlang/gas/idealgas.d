/**
 * idealgas.d
 * Ideal gas model for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22: initial cut, to explore options.
 */

module idealgas;
import gasmodel;
import perfectgasEOS;
import calperfectgasEOS;
import sutherland_visc;
import sutherland_therm_cond;
import std.math;
import std.stdio;
import std.file;
import std.json;
import std.conv;
import luad.all;
import std.c.stdlib : exit;

class IdealGas: GasModel {
public:
    this() {
	// Default model is mostly initialized in the private data below.
	_n_species = 1;
	_n_modes = 1;
	_species_names ~= "ideal air";
    }
    this(in string species_name, in double[string] params) {
	this();
	// The species name was set to 'ideal air' by default.
	// Let's overwrite that here.
	_species_names[0] = species_name;
	// Now, pull out the remaining numeric value parameters.
	_Mmass = params["Mmass"];
	_gamma = params["gamma"];
	// Reference values for entropy
	_s1 = params["s1"];
	_T1 = params["T1"];
	_p1 = params["p1"];
	// Molecular transport coefficent constants.
	_mu_ref = params["mu_ref"];
	_T_ref = params["T_ref"];
	_S_mu = params["S_mu"];
	_k_ref = params["k_ref"];
	_S_k = params["S_k"];
	// Compute derived parameters
	_Rgas = R_universal/_Mmass;
	_Cv = _Rgas / (_gamma - 1.0);
	_Cp = _Rgas*_gamma/(_gamma - 1.0);
    }

    override string toString() const
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

    override void update_thermo_from_pT(ref GasState Q) const 
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.rho = density(Q.p, Q.T[0], _Rgas);
	Q.e[0] = energy(Q.T[0], _Cv, 0.0);
    }
    override void update_thermo_from_rhoe(ref GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.T[0] = calperfectgasEOS.temperature(Q.e[0], _Cv, 0.0);
	Q.p = pressure(Q.rho, Q.T[0], _Rgas);
    }
    override void update_thermo_from_rhoT(ref GasState Q) const
    {
	throw new Exception("not implemented");
    }
    override void update_thermo_from_rhop(ref GasState Q) const
    {
	throw new Exception("not implemented");
    }
    override void update_thermo_from_ps(ref GasState Q, double s) const
    {
	throw new Exception("not implemented");
    }
    override void update_thermo_from_hs(ref GasState Q, double s) const
    {
	throw new Exception("not implemented");
    }
    override void update_sound_speed(ref GasState Q) const
    {
	Q.a = sqrt(_gamma * _Rgas * Q.T[0]);
    }
    override void update_trans_coeffs(ref GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect number of modes");
	assert(Q.k.length == 1, "incorrect number of modes");
	Q.mu = sutherland_viscosity(Q.T[0], _T_ref, _mu_ref, _S_mu);
	Q.k[0] = sutherland_thermal_conductivity(Q.T[0], _T_ref, _k_ref, _S_k);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
	throw new Exception("not implemented");
    }
    */
    override double dedT_const_v(in GasState Q) const
    {
	return _Cv;
    }
    override double dhdT_const_p(in GasState Q) const
    {
	return _Cp;
    }
    override double gas_constant(in GasState Q) const
    {
	return R_universal/_Mmass;
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.e[0];
    }
    override double enthalpy(in GasState Q) const
    {
	return Q.e[0] + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
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

} // end class Ideal_gas

IdealGas init_ideal_gas(ref LuaTable t)
{
    auto species_name = t.get!string("species_name");
    string[10] plist = ["Mmass", "gamma",
			"s1", "T1", "p1",
			"mu_ref", "T_ref", "S_mu", "k_ref", "S_k"];
    double[string] params;
    foreach (p; plist) {
	try {
	    params[p] = t.get!double(p);
	} catch (Exception e) {
	    writeln("ERROR: There was a problem reading the value for ", p);
	    writeln("ERROR: when initialising the ideal gas model.");
	    writeln("ERROR: Quitting at this point.");
	    exit(1);
	}
    }

    auto gm = new IdealGas(species_name, params);
    return gm;
}

unittest {
    import std.stdio;
    auto gm = new IdealGas();
    assert(gm.species_name(0) == "ideal air", "species name list");
    auto gd = GasState(gm, 100.0e3, 300.0);
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

    auto lua = new LuaState;
    lua.openLibs();
    lua.doFile("sample-data/ideal-air-gas-model.lua");
    auto t = lua.get!LuaTable("ideal_gas");
    gm = init_ideal_gas(t);
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
