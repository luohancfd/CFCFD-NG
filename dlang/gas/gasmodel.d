/**
 * gasmodel.d
 * Storage arrangement for the data defining a gas state,
 * interface description of the gas model functionality and
 * utilities to create specific gas model objects.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22, first cut, exploring the options.
 */

module gasmodel;
import std.conv;

immutable double R_universal = 8.314; // J/mole.K

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
    const void update_thermo_from_hs(ref GasState Q, double s) {}
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


class GasState {
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

    this() {} // makes no sense to define the data in the absence of a model

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

    // Note that we must not send the current object in the others list as well.
    void copy_average_values_from(in GasState[] others, in GasModel gm) 
    {
	uint n = others.length;
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
    }

    override string toString()
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
