/**
 * gas_model.d
 * Storage arrangement for the data defining a gas state,
 * interface description of the gas model functionality and
 * utilities to create specific gas model objects.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22, first cut, exploring the options.
 */

module gas_model;

immutable double R_universal = 8.314; // J/mole.K

class Gas_model {
public:
    @property uint n_species() { return _n_species; }
    @property uint n_modes() { return _n_modes; }

    // Methods to be overridden.
    void eval_thermo_state_pT(ref Gas_data Q) {}
    void eval_thermo_state_rhoe(ref Gas_data Q) {}
    void eval_thermo_state_rhoT(ref Gas_data Q) {}
    void eval_thermo_state_rhop(ref Gas_data Q) {}
    void eval_thermo_state_ps(ref Gas_data Q, double s) {}
    void eval_thermo_state_hs(ref Gas_data Q, double s) {}
    void eval_sound_speed(ref Gas_data Q) {}
    void eval_transport_coefficients(ref Gas_data Q) {}
    void eval_diffusion_coefficients(ref Gas_data Q) {}

    // Methods to be overridden.
    double dedT_const_v(in Gas_data Q) { return 0.0; }
    double dhdT_const_p(in Gas_data Q) { return 0.0; }
    double gas_constant(in Gas_data Q) { return 0.0; }
    double internal_energy(in Gas_data Q) { return 0.0; }
    double enthalpy(in Gas_data Q) { return 0.0; }
    double entropy(in Gas_data Q) { return 0.0; }
    
    final double Cv(in Gas_data Q) { return dedT_const_v(Q); }
    final double Cp(in Gas_data Q) { return dhdT_const_p(Q); }
    final double R(in Gas_data Q) { return gas_constant(Q); }
    final double gamma(in Gas_data Q) { return Cp(Q)/Cv(Q); }
protected:
    uint _n_species;
    uint _n_modes;
}


class Gas_data {
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
    double[][] D_AB; /// binary diffusion coefficients
    double sigma;    /// electrical conductivity, S/m
    /// Composition
    double[] massf;  /// species mass fractions
    double quality;  /// vapour quality

    this(Gas_model gm, double p_init=100.0e3, double T_init=300.0, 
	 double[] massf_init=[1.0,], double quality_init=1.0)
    {
	p = p_init;
	p_e = p_init;
	T.length = gm.n_modes;
	foreach(ref Tmode; T) {
	    Tmode = T_init;
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
	// Now, evaluate the rest of the properties using the gas model.
	gm.eval_thermo_state_pT(this);
	gm.eval_sound_speed(this);
	gm.eval_transport_coefficients(this);
    }

    this() {} // makes no sense to define the data in the absence of a model

    this(Gas_data other) {
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
} // end class Gas_data
