/**
 * gasmodel.d
 * Storage arrangement for the data defining a gas state,
 * interface description of the gas model functionality and
 * utilities to create specific gas model objects.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22, first cut, exploring the options.
 */

module gas.gas_model;

import std.conv;
import std.math;
import std.stdio;
import std.range;
import std.algorithm;

import gas.physical_constants;

immutable double SMALL_MOLE_FRACTION = 1.0e-15;
immutable double T_MIN = 20.0; 

class GasModel {
public:
    @property const uint n_species() { return _n_species; }
    @property const uint n_modes() { return _n_modes; }
    const string species_name(int i) { return _species_names[i]; }

    // Methods to be overridden.
    const void update_thermo_from_pT(ref GasState Q) {}
    const void update_thermo_from_rhoe(ref GasState Q) {}
    const void update_thermo_from_rhoT(ref GasState Q) {}
    const void update_thermo_from_rhop(ref GasState Q) {}
    const void update_thermo_from_ps(ref GasState Q, double s) {}
    const void update_thermo_from_hs(ref GasState Q, double h, double s) {}
    const void update_sound_speed(ref GasState Q) {}
    const void update_trans_coeffs(ref GasState Q) {}
    // const void update_diff_coeffs(ref GasState Q) {}

    // Methods to be overridden.
    const double dedT_const_v(in GasState Q) { return 0.0; }
    const double dhdT_const_p(in GasState Q) { return 0.0; }
    const double gas_constant(in GasState Q) { return 0.0; }
    const double internal_energy(in GasState Q) { return 0.0; }
    const double enthalpy(in GasState Q) { return 0.0; }
    const double entropy(in GasState Q) { return 0.0; }
    
    final const double Cv(in GasState Q) { return dedT_const_v(Q); }
    final const double Cp(in GasState Q) { return dhdT_const_p(Q); }
    final const double R(in GasState Q) { return gas_constant(Q); }
    final const double gamma(in GasState Q) { return Cp(Q)/Cv(Q); }
protected:
    // These data need to be properly initialized by the derived class.
    uint _n_species;
    uint _n_modes;
    string[] _species_names;
}


struct GasState {
public:
    /// Thermodynamic properties.
    double rho;  /// density, kg/m**3
    double p;    /// presure, Pa
    double p_e;  /// electron pressure
    double a;    /// sound speed, m/s
    double[] e;  /// specific internal energies, J/kg
    double[] T;  /// temperatures, K
    /// Transport properties
    double mu;   /// viscosity, Pa.s
    double[] k;  /// thermal conductivities, W/(m.k)
    // double[][] D_AB; /// binary diffusion coefficients
    double sigma;    /// electrical conductivity, S/m
    /// Composition
    double[] massf;  /// species mass fractions
    double quality;  /// vapour quality

    this(uint n_species, uint n_modes)
    {
	massf.length = n_species;
	e.length = n_modes;
	T.length = n_modes;
	k.length = n_modes;
    }

    this(in GasModel gm, in double p_init, in double[] T_init, 
	 in double[] massf_init=[1.0,], in double quality_init=1.0,
	 in double sigma_init=0.0)
    {
	p = p_init;
	p_e = p_init;
	T.length = gm.n_modes;
	foreach(i; 0 .. gm.n_modes) {
	    try
		{ T[i] = T_init[i]; }
	    catch (Exception e)
		// We assume that at least 1 temperature arrived.
		{ T[i] = T_init[0]; }
	}
	e.length = gm.n_modes;
	k.length = gm.n_modes;
	massf.length = gm.n_species;
	foreach(i; 0 .. gm.n_species) {
	    try 
		{ massf[i] = massf_init[i]; }
	    catch (Exception e) 
		{ massf[i] = 0.0; }
	}
	quality = quality_init;
	sigma = sigma_init;
	// Now, evaluate the rest of the properties using the gas model.
	gm.update_thermo_from_pT(this);
	gm.update_sound_speed(this);
	gm.update_trans_coeffs(this);
    }

    this(in GasModel gm, in double p_init, in double T_init, 
	 in double[] massf_init=[1.0,], in double quality_init=1.0,
	 in double sigma_init=0.0)
    {
	double[] Tlocal;
	Tlocal.length = gm.n_modes;
	foreach(ref Tmode; Tlocal) {
	    Tmode = T_init;
	}
	this(gm, p_init, Tlocal, massf_init, quality_init, sigma_init);
    }

    this(in GasState other) 
    {
	rho = other.rho;
	p = other.p;
	p_e = other.p_e;
	a = other.a;
	e = other.e.dup;
	T = other.T.dup;
	mu = other.mu;
	k = other.k.dup;
	// D_AB
	sigma = other.sigma;
	massf = other.massf.dup;
	quality = other.quality;
    }

    // Postblit constructor
    this(this)
    {
	massf = massf.dup;
	e = e.dup;
	T = T.dup;
	k = k.dup;
    }

    void copy_values_from(in GasState other) 
    {
	rho = other.rho;
	p = other.p;
	p_e = other.p_e;
	a = other.a;
	e = other.e.dup;
	T = other.T.dup;
	mu = other.mu;
	k = other.k.dup;
	// D_AB
	sigma = other.sigma;
	massf = other.massf.dup;
	quality = other.quality;
    }

    void copy_average_values_from(in GasState gs0, in GasState gs1) 
    // Avoids memory allocation, it's all in place.
    {
	rho = 0.5 * (gs0.rho + gs1.rho);
	p = 0.5 * (gs0.p + gs1.p);
	p_e = 0.5 * (gs0.p_e + gs1.p_e);
	a = 0.5 * (gs0.a + gs1.a);
	foreach(i; 0 .. e.length) e[i] = 0.5 * (gs0.e[i] + gs1.e[i]);
	foreach(i; 0 .. T.length) T[i] = 0.5 * (gs0.T[i] + gs1.T[i]);
	mu = 0.5 * (gs0.mu + gs1.mu);
	foreach(i; 0 .. k.length) k[i] = 0.5 * (gs0.k[i] + gs1.k[i]);
	// D_AB
	sigma = 0.5 * (gs0.sigma * gs1.sigma);
	foreach(i; 0 .. massf.length) massf[i] = 0.5 * (gs0.massf[i] * gs1.massf[i]);
	quality = 0.5 * (gs0.quality + gs1.quality);
    }

    void copy_average_values_from(in GasState[] others, in GasModel gm) 
    // Note that we must not send the current object in the others list as well.
    {
	size_t n = others.length;
	if (n == 0) throw new Error("Need to average from a nonempty array.");
	foreach(other; others) {
	    if ( this is other ) throw new Error("Must not include destination in source list.");
	}
	// Accumulate from a clean slate and then divide.
	p = 0.0;
	p_e = 0.0;
	foreach(ref elem; T) elem = 0.0;
	sigma = 0.0;
	foreach(ref elem; massf) elem = 0.0;
	quality = 0.0;
	foreach(other; others) {
	    p += other.p;
	    p_e += other.p_e;
	    foreach(i; 0 .. T.length) T[i] += other.T[i];
	    sigma += other.sigma;
	    foreach(i; 0 .. massf.length) massf[i] += other.massf[i];
	    quality += other.quality;
	}
	p /= n;
	p_e /= n;
	foreach(ref elem; T) elem /= n;
	sigma /= n;
	foreach(ref elem; massf) elem /= n;
	quality /= n;
	// Now, evaluate the rest of the properties using the gas model.
	gm.update_thermo_from_pT(this);
	gm.update_sound_speed(this);
	gm.update_trans_coeffs(this);
    } // end copy_average_values_from()

    bool check_values(bool print_message=true) const
    {
	double RHOMIN = 0.0;
	double TMIN = 1.0;
	bool is_data_valid = true;

	if ( !(isFinite(rho)) || rho < 1.01 * RHOMIN ) {
	    if (print_message) writeln("Density invalid: ", rho);
	    is_data_valid = false;
	}
	auto nmodes = e.length;
	foreach(imode; 0 .. nmodes) {
	    if ( !isFinite(T[imode]) || T[imode] < 1.01 * TMIN ) {
		if ( print_message ) writeln("Temperature[", imode, "] invalid: ", T[imode]);
		is_data_valid = false;
	    }
	    if ( !isFinite(e[imode]) ) {
		if ( print_message ) writeln("Energy[", imode, "] invalid: ", e[imode]);
		is_data_valid = false;
	    }
	}
	if ( !isFinite(p) ) {
	    if ( print_message ) writeln("Total pressure invalid: ", p);
	    is_data_valid = false;
	}
	if ( !isFinite(p_e) ) {
	    if ( print_message ) writeln("Electron pressure invalid: ", p_e);
	    is_data_valid = false;
	}
	if ( !isFinite(a) ) {
	    if ( print_message ) writeln("Sound speed invalid: ", a);
	    is_data_valid = false;
	}
	double f_sum = 0.0; foreach(elem; massf) f_sum += elem;
	if ( f_sum < 0.99 || f_sum > 1.01 || !isFinite(f_sum) ) {
	    if ( print_message ) writeln("Mass fraction sum bad: ", f_sum);
	    is_data_valid = false;
	}
	return is_data_valid;
    } // end check_values()

    string toString() const
    {
	char[] repr;
	repr ~= "GasState(";
	repr ~= "rho=" ~ to!string(rho);
	repr ~= ", p=" ~ to!string(p);
	repr ~= ", p_e=" ~ to!string(p_e);
	repr ~= ", a=" ~ to!string(a);
	repr ~= ", T=" ~ to!string(T);
	repr ~= ", e=" ~ to!string(e);
	repr ~= ", mu=" ~ to!string(mu);
	repr ~= ", k=" ~ to!string(k);
	repr ~= ", massf=" ~ to!string(massf);
	repr ~= ", quality=" ~ to!string(quality);
	repr ~= ", sigma=" ~ to!string(sigma);
	repr ~= ")";
	return to!string(repr);
    }

/+ [TODO]
    double * copy_values_to_buffer(double *buf) const;
    double * copy_values_from_buffer(double *buf);
+/
} // end class GasState


void scale_mass_fractions(ref double[] massf, double tolerance=0.0)
{
    auto my_nsp = massf.length;
    if ( my_nsp == 1 ) {
	// Single species, always expect massf[0]==1.0, so we can take a short-cut.
	if ( fabs(massf[0] - 1.0) > 0.1 ) 
	    throw new Error("Single species mass fraction far from 1.0");
	massf[0] = 1.0;
    } else {
	// Multiple species, do the full job.
	double massf_sum = 0.0;
	foreach(isp; 0 .. my_nsp) {
	    massf[isp] = massf[isp] >= 0.0 ? massf[isp] : 0.0;
	    massf_sum += massf[isp];
	}
	if ( fabs(massf_sum - 1.0) > 0.1 )
	    throw new Error("Sum of species mass fractions far from 1.0");
	if ( fabs(massf_sum - 1.0) > tolerance ) {
	    foreach(isp; 0 .. my_nsp) massf[isp] /= massf_sum;
	}
    }
} // end scale_mass_fractions()

pure double mass_average(in GasState Q, in double[] phi)
{
    assert(Q.massf.length == phi.length);
    return reduce!((acc, e) =>  acc + e[0]*e[1])(0.0, zip(Q.massf, phi));
}

double mixture_molecular_weight(in double[] massf, in double[] MW)
{
    assert(massf.length == MW.length);
    double M_inv = reduce!((acc, e) => acc + e[0]/e[1])(0.0, zip(massf, MW));
    return 1.0/M_inv;
}

void massf2molef(in double[] massf, in double[] MW, double[] molef)
{
    assert(massf.length == MW.length);
    assert(massf.length == molef.length);
    double MW_mix = mixture_molecular_weight(massf, MW);
    for ( auto isp = 0; isp < massf.length; ++isp ) {
	molef[isp] = massf[isp]*MW_mix/MW[isp];
    }
}

unittest {
    auto gd = GasState(2, 1);
    gd.massf[0] = 0.8;
    gd.massf[1] = 0.2;
    double[] phi = [9.0, 16.0];
    assert(approxEqual(10.4, mass_average(gd, phi)));
    double[] MW = [28.0e-3, 32.0e-3];
    assert(approxEqual(0.0287179, mixture_molecular_weight(gd.massf, MW)));
    double[] molef = [0.0, 0.0];
    massf2molef(gd.massf, MW, molef);
    assert(approxEqual([0.8205128, 0.1794872], molef));
}
