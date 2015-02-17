/**
 * therm_perf_gas_mix.d
 * Thermally perfect gas mix model for use in the CFD codes.
 *
 * Author: Rowan G. and Peter J.
 * First code: 27-Jan-2015
 */

module gas.therm_perf_gas;

import gas.gas_model;
import gas.physical_constants;
import gas.thermo.cea_thermo_curves;
import gas.thermo.perf_gas_mix_eos;
import gas.thermo.therm_perf_gas_mix_eos;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.diffusion.cea_viscosity;
import gas.diffusion.cea_therm_cond;
import gas.diffusion.wilke_mixing_viscosity;
import gas.diffusion.wilke_mixing_therm_cond;
import std.math;
import std.stdio;
import std.string;
import luad.all;
import util.lua_service;
import ridder;
import std.c.stdlib : exit;

class ThermallyPerfectGas: GasModel {
public:
    this(in string fname)
    {
	auto lua = initLuaState(fname);
	getArray!string(lua.get!LuaTable("species"), _species_names, "species");
	_n_species = cast(uint) _species_names.length;
	_n_modes = 1;
	
	// 1. Initialise gas constants from molecular mass
	_R.length = _n_species;
	double[] Mmass;
	Mmass.length = _n_species;
	foreach ( isp; 0.._n_species ) {
	    Mmass[isp] = lua.get!double(_species_names[isp], "M");
	    _R[isp] = R_universal/Mmass[isp];
	}
	// 2. Set the thermal EOS
	_thermalEOS = new PerfectGasMixEOS(_R);
	// 3. Set the caloric EOS
	foreach ( isp; 0.._n_species ) {
	    _curves ~= createCEAThermo(lua.get!LuaTable(_species_names[isp], "cea_thermo"), _R[isp]);
	}
	_caloricEOS = new ThermallyPerfectGasMixEOS(_R, _curves);
	// 4. Set the viscosity model
	Viscosity[] vms;
	foreach ( isp; 0.._n_species ) {
	    vms ~= createCEAViscosity(lua.get!LuaTable(_species_names[isp], "viscosity"));
	}
	_viscModel = new WilkeMixingViscosity(vms, Mmass);
	// 5. Set the thermal conductivity model
	ThermalConductivity[] tcms;
	foreach (isp; 0.._n_species ) {
	    tcms ~= createCEAThermalConductivity(lua.get!LuaTable(_species_names[isp], "therm_cond"));

	}
	_thermCondModel = new WilkeMixingThermCond(tcms, Mmass);
    }

    override string toString() const
    {
	return "";
    }

    @nogc override void update_thermo_from_pT(ref GasState Q) const 
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	_thermalEOS.update_density(Q);
	_caloricEOS.update_energy(Q);
    }
    @nogc override void update_thermo_from_rhoe(ref GasState Q) const
    {
	assert(Q.e.length == 1, "incorrect length of energy array");
	_caloricEOS.update_temperature(Q);
	_thermalEOS.update_pressure(Q);
    }
    @nogc override void update_thermo_from_rhoT(ref GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	_caloricEOS.update_energy(Q);
	_thermalEOS.update_pressure(Q);
    }
    @nogc override void update_thermo_from_rhop(ref GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	_thermalEOS.update_temperature(Q);
	_caloricEOS.update_energy(Q);
    }
    override void update_thermo_from_ps(ref GasState Q, double s) const
    {
	double TOL = 1.0e-6;
	double delT = 100.0;
	double T1 = Q.T[0];
	double Tsave = T1;
	double T2 = T1 + delT;

	auto zeroFun = delegate (double T) {
	    Q.T[0] = T;
	    double s_guess = entropy(Q);
	    return s - s_guess;
	};

	if ( bracket!zeroFun(T1, T2) == -1 ) {
	    string msg = "The 'bracket' function failed to find temperature values\n";
	    msg ~= "that bracketed the zero function in ThermallyPerfectGas.update_thermo_from_ps().\n";
	    msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
	    throw new Exception(msg);
	}

	if ( T1 < T_MIN )
	    T1 = T_MIN;
	
	try {
	    Q.T[0] = solve!zeroFun(T1, T2, TOL);
	}
	catch ( Exception e ) {
	    string msg = "There was a problem iterating to find temperature\n";
	    msg ~= "in function ThermallyPerfectGas.update_thermo_from_ps().\n";
	    msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
	    msg ~= format("The target entropy value was: %12.6f\n", s);
	    msg ~= format("The GasState is currently:\n");
	    msg ~= Q.toString();
	    msg ~= "The message from the ridder.solve function is:\n";
	    msg ~= e.msg;
	    throw new Exception(msg);
	}
	_caloricEOS.update_energy(Q);
	_thermalEOS.update_density(Q);
    }
    override void update_thermo_from_hs(ref GasState Q, double h, double s) const
    {
	// We do this in two stages.
	// First, from enthalpy we compute temperature.
	double TOL = 1.0e-6;
	double delT = 100.0;
	double T1 = Q.T[0];
	double Tsave = T1;
	double T2 = T1 + delT;

	auto zeroFun = delegate (double T) {
	    Q.T[0] = T;
	    double h_guess = enthalpy(Q);
	    return h - h_guess;
	};

	if ( bracket!zeroFun(T1, T2) == -1 ) {
	    string msg = "The 'bracket' function failed to find temperature values\n";
	    msg ~= "that bracketed the zero function in ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
	    throw new Exception(msg);
	}

	if ( T1 < T_MIN )
	    T1 = T_MIN;
	
	try {
	    Q.T[0] = solve!zeroFun(T1, T2, TOL);
	}
	catch ( Exception e ) {
	    string msg = "There was a problem iterating to find temperature\n";
	    msg ~= "in function ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
	    msg ~= format("The target enthalpy value was: %12.6f\n", h);
	    msg ~= format("The GasState is currently:\n");
	    msg ~= Q.toString();
	    msg ~= "The message from the ridder.solve function is:\n";
	    msg ~= e.msg;
	    throw new Exception(msg);
	}

	// Second, we can iterate to find the pressure that gives
	// correct entropy.
	TOL = 1.0e-3;
	double delp = 1000.0;
	double p1 = Q.p;
	double psave = p1;
	double p2 = p1 + delp;

	auto zeroFun2 = delegate (double p) {
	    Q.p = p;
	    double s_guess = entropy(Q);
	    return s - s_guess;
	};

	if ( bracket!zeroFun2(p1, p2) == -1 ) {
	    string msg = "The 'bracket' function failed to find pressure values\n";
	    msg ~= "that bracketed the zero function in ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The final values are: p1 = %12.6f and p2 = %12.6f\n", p1, p2);
	    throw new Exception(msg);
	}

	if ( p1 < 0.0 )
	    p1 = 1.0;
	
	try {
	    Q.p = solve!zeroFun2(p1, p2, TOL);
	}
	catch ( Exception e ) {
	    string msg = "There was a problem iterating to find pressure\n";
	    msg ~= "in function ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The initial pressure guess was: %12.6f\n", psave);
	    msg ~= format("The target entropy value was: %12.6f\n", s);
	    msg ~= format("The GasState is currently:\n");
	    msg ~= Q.toString();
	    msg ~= "The message from the ridder.solve function is:\n";
	    msg ~= e.msg;
	    throw new Exception(msg);
	}
	_caloricEOS.update_energy(Q);
	_thermalEOS.update_density(Q);
    }
    @nogc override void update_sound_speed(ref GasState Q) const
    {
	Q.a = double.init;
    }
    override void update_trans_coeffs(ref GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect number of modes");
	assert(Q.k.length == 1, "incorrect number of modes");
	_viscModel.update_viscosity(Q);
	_thermCondModel.update_thermal_conductivity(Q);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
	throw new Exception("not implemented");
    }
    */
    override double dedT_const_v(in GasState Q) const
    {
	throw new Exception("not implemented");
    }
    override double dhdT_const_p(in GasState Q) const
    {
	throw new Exception("not implemented");
    }
    override double gas_constant(in GasState Q) const
    {
	return mass_average(Q, _R);
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.e[0];
    }
    override double enthalpy(in GasState Q) const
    {
	double[] h = new double[_n_species];
	foreach ( isp; 0.._n_species) {
	    h[isp] = _curves[isp].eval_h(Q.T[0]);
	}
	return mass_average(Q, h);
    }
    override double entropy(in GasState Q) const
    {
	double[] s = new double[_n_species];
	foreach ( isp; 0.._n_species ) {
	    s[isp] = _curves[isp].eval_s(Q.T[0]) - _R[isp]*log(Q.p/P_atm);
	}
	return mass_average(Q, s);
    }

private:
    double[] _R;
    PerfectGasMixEOS _thermalEOS;
    ThermallyPerfectGasMixEOS _caloricEOS;
    CEAThermo[] _curves;
    WilkeMixingViscosity _viscModel;
    WilkeMixingThermCond _thermCondModel;
} // end class Ideal_gas

unittest 
{
    import util.msg_service;

    auto gm = new ThermallyPerfectGas("sample-data/therm-perf-5-species-air.lua");
    auto gd = new GasState(5, 1);
    gd.p = 1.0e6;
    gd.T[0] = 2000.0;
    gd.massf = [0.2, 0.2, 0.2, 0.2, 0.2];
    gm.update_thermo_from_pT(gd);
    assert(approxEqual(11801825.6, gd.e[0]), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(1.2840117, gd.rho), failedUnitTest(__LINE__, __FILE__));

    gd.rho = 2.0;
    gd.e[0] = 14.0e6;
    gm.update_thermo_from_rhoe(gd);
    assert(approxEqual(3373757.4, gd.p), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(4331.944, gd.T[0]), failedUnitTest(__LINE__, __FILE__));
    
    gd.T[0] = 10000.0;
    gd.rho = 1.5;
    gm.update_thermo_from_rhoT(gd);
    assert(approxEqual(5841068.3, gd.p), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(20340105.9, gd.e[0]), failedUnitTest(__LINE__, __FILE__));

    gd.rho = 10.0;
    gd.p = 5.0e6;
    gm.update_thermo_from_rhop(gd);
    assert(approxEqual(11164648.5, gd.e[0]), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(1284.012, gd.T[0]), failedUnitTest(__LINE__, __FILE__));

    gd.p = 1.0e6;
    double s = 10000.0;
    gm.update_thermo_from_ps(gd, s);
    assert(approxEqual(2560.118, gd.T[0]), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(12313952.52, gd.e[0]), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(1.00309, gd.rho), failedUnitTest(__LINE__, __FILE__));

    s = 11000.0;
    double h = 17.0e6;
    gm.update_thermo_from_hs(gd, h, s);
    assert(approxEqual(5273.103, gd.T[0]), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(14946629.7, gd.e[0]), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(0.4603513, gd.rho), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(945271.84, gd.p), failedUnitTest(__LINE__, __FILE__));

    gd.T[0] = 4000.0;
    gm.update_trans_coeffs(gd);
    assert(approxEqual(0.00012591, gd.mu), failedUnitTest(__LINE__, __FILE__));
    assert(approxEqual(0.2448263, gd.k[0]), failedUnitTest(__LINE__, __FILE__));

    // TODO: entropy, enthalpy tests.
}
