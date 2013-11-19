/** \file gas-model.hh
 *  \ingroup gas
 *
 *  \author Rowan J Gollan
 *  \version 03-Jul-2008
 **/

#ifndef GAS_MODEL_HH
#define GAS_MODEL_HH

#include <iostream>
#include <vector>
#include <map>
#include <valarray>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../nm/source/functor.hh"
#include "gas_data.hh"

// Some constants for the module
#define DEFAULT_MIN_MOLE_FRACTION 1.0e-15
#define DEFAULT_MIN_MASS_FRACTION 1.0e-10

const double rho_min = 1.0e-10;
const double T_min = 1.0e-9;
const double T_max = 150000.0;
const double p_min = 1.0e-6;
const double a_min = 1.0e-2;

/**
 * \brief Class for encapsulating gas model behaviour.
 *
 * The Gas_model class is an abstract base class. It 
 * prescribes the services offered by all derived classes.
 * This class encapsulates the behaviour of the gas model.
 * Basically, given a Gas_data object, a Gas_model object
 * can compute various properties based on the model behaviour.
 * It is a catch-all class which presents a single, unified
 * object for the calling code to access all gas property-related
 * computations. These computations include:
 *  - evaluation of thermodynamic state and derivatives
 *  - conversion to/from conserved quantities and primary variables
 *  - evaluation of transport properties
 *  - evaluation of general gas properties
 *  - evaluation of gas properties of specific components
 **/

class Gas_model {
public:
    Gas_model();
    Gas_model(std::string cfile);
    virtual ~Gas_model();

    int get_number_of_species() const
    { return nsp_; }
    
    int get_number_of_modes() const
    { return nmodes_; }

    void set_reaction_compatibility(bool val)
    { good_for_reactions_ = val; }

    bool good_for_reactions() const
    { return good_for_reactions_; }

    int number_of_values_in_gas_data_copy() const;

    double gas_rhomin()
    { return rho_min; }
    
    double gas_Tmin()
    { return T_min; }

    double gas_Tmax()
    { return T_max; }

    double gas_pmin()
    { return p_min; }
    
    double gas_amin()
    { return a_min; }

    int eval_thermo_state_pT(Gas_data &Q)
    { return s_eval_thermo_state_pT(Q); }

    int eval_thermo_state_rhoe(Gas_data &Q) 
    { return s_eval_thermo_state_rhoe(Q); }

    int eval_thermo_state_rhoT(Gas_data &Q)
    { return s_eval_thermo_state_rhoT(Q); }

    int eval_thermo_state_rhop(Gas_data &Q)
    { return s_eval_thermo_state_rhop(Q); }

    int eval_thermo_state_ps(Gas_data &Q, double p, double s)
    { return s_eval_thermo_state_ps(Q, p, s); }

    int eval_thermo_state_hs(Gas_data &Q, double h, double s)
    { return s_eval_thermo_state_hs(Q, h, s); }

    int eval_sound_speed(Gas_data &Q)
    { return s_eval_sound_speed(Q); }

    int eval_transport_coefficients(Gas_data &Q)
    { return s_eval_transport_coefficients(Q); }

    int eval_diffusion_coefficients(Gas_data &Q)
    { return s_eval_diffusion_coefficients(Q); }

    // State derivatives

    double dTdp_const_rho(const Gas_data &Q, int &status)
    { return s_dTdp_const_rho(Q, status); }
    
    double dTdp_const_rho_py(const Gas_data &Q)
    { int status; return s_dTdp_const_rho(Q, status); }

    double dTdrho_const_p(const Gas_data &Q, int &status)
    { return s_dTdrho_const_p(Q, status); }

    double dTdrho_const_p_py(const Gas_data &Q)
    { int status; return s_dTdrho_const_p(Q, status); }
    
    double dTdrho_const_s_py( const Gas_data &Q )
    { int status; return s_dTdrho_const_s( Q, status ); }

    double dpdrho_const_T(const Gas_data &Q, int &status)
    { return s_dpdrho_const_T(Q, status); }
    
    double dpdrho_i_const_T(const Gas_data &Q, int isp, int &status)
    { return s_dpdrho_i_const_T(Q, isp, status); }

    double dedT_const_v(const Gas_data &Q, int &status)
    { return s_dedT_const_v(Q, status); }
    
    double dhdT_const_p(const Gas_data &Q, int &status)
    { return s_dhdT_const_p(Q, status); }
    
    double dpdT_i_const_rho(const Gas_data &Q, int itm, int &status)
    { return s_dpdT_i_const_rho(Q, itm, status); }
    
    // Properties of the gas
    // alias for dedT_const_v
    double Cv(const Gas_data &Q, int &status)
    { return s_dedT_const_v(Q, status); }

    double Cp(const Gas_data &Q, int &status)
    { return s_dhdT_const_p(Q, status); }

    double R(const Gas_data &Q, int &status)
    { return s_gas_constant(Q, status); }

    double gamma(const Gas_data &Q, int &status)
    { return Cp(Q, status)/Cv(Q, status); }

    // Computed from definition
    double Prandtl(double mu, double Cp, double k)
    { return mu*Cp/k; }

    double mixture_molecular_weight(const Gas_data &Q)
    { return s_mixture_molecular_weight(Q); }

    // Some component properties
    double molecular_weight(int isp)
    { return s_molecular_weight(isp); }

    double internal_energy(const Gas_data &Q, int isp)
    { return s_internal_energy(Q, isp); }

    double enthalpy(const Gas_data &Q, int isp)
    { return s_enthalpy(Q, isp); }
    
    double modal_enthalpy( const Gas_data &Q, int isp, int itm)
    { return s_modal_enthalpy(Q,isp,itm); }
    
    double modal_Cv( Gas_data &Q, int itm )
    { return s_modal_Cv(Q,itm); }

    double entropy(const Gas_data &Q, int isp)
    { return s_entropy(Q, isp); }

    double Gibbs_free_energy(const Gas_data &Q, int isp);

    std::string species_name(int isp)
    { return s_names_[isp]; }
    
    int get_isp_from_species_name( std::string name );

    std::string mode_name(int imode)
    { return m_names_[imode]; }

    std::string mode_component_name(int imode, int ic)
    { return m_components_[imode][ic]; }

    int mode_no_components(int imode)
    { return m_components_[imode].size(); }

    int get_imode_from_mode_name(std::string name);
    
    int set_mole_fractions(Gas_data &Q, 
			   std::vector<std::string> &sp, 
			   const std::vector<double> &X);

    int no_atoms_of(std::string atom, int isp);
    void atomic_constituents(int isp, std::map<std::string, int> &m);
    int charge(int isp)
    { return charge_[isp]; }
    
    // Some useful analysis functions
    double mixture_internal_energy(const Gas_data &Q, double T=0.0);
    double mixture_enthalpy(const Gas_data &Q, double T=0.0);
    double mixture_entropy(const Gas_data &Q);
    
    // Overloaded functions for python
    double R(const Gas_data &Q)
    { int status; return s_gas_constant(Q, status); }
    
    double Cv(const Gas_data &Q)
    { int status; return s_dedT_const_v(Q, status); }
    
    double Cp(const Gas_data &Q)
    { int status; return s_dhdT_const_p(Q, status); }
    
    double gamma(const Gas_data &Q)
    { int status; return Cp(Q, status)/Cv(Q, status); }

    const std::vector<double>& M() { return M_; }


protected:
    int nsp_;     // No. of species components
    int nmodes_;  // No. of (separate) thermal modes
    bool good_for_reactions_; // indicates if the model can be used with the
                              // finite-rate chemistry module

    std::vector<std::string> s_names_;
    std::vector<std::string> m_names_;
    std::vector<std::vector<std::string> > m_components_;
    std::vector<double> M_;
    std::vector<int> charge_;
    std::vector<std::map<std::string, int> > atomic_constituents_;

    // Derived clasess may set nsp_ and nmodes_
    void set_number_of_species(int nsp)
    { nsp_ = nsp; }

    void set_number_of_modes(int nmodes)
    { nmodes_ = nmodes; }

    // Derived classes need to implement their own versions of the following...
    virtual int s_eval_thermo_state_rhoe(Gas_data &Q) = 0;
    virtual int s_eval_thermo_state_pT(Gas_data &Q);
    virtual int s_eval_thermo_state_rhoT(Gas_data &Q);
    virtual int s_eval_thermo_state_rhop(Gas_data &Q);
    virtual int s_eval_thermo_state_ps(Gas_data &Q, double p, double s);
    virtual int s_eval_thermo_state_hs(Gas_data &Q, double h, double s);
    virtual int s_eval_sound_speed(Gas_data &Q);
    virtual int s_eval_transport_coefficients(Gas_data &Q) = 0;
    virtual int s_eval_diffusion_coefficients(Gas_data &Q) = 0;
    virtual double s_dTdp_const_rho(const Gas_data &Q, int &status);
    virtual double s_dTdrho_const_p(const Gas_data &Q, int &status);
    virtual double s_dTdrho_const_s(const Gas_data &Q, int &status );
    virtual double s_dpdrho_const_T(const Gas_data &Q, int &status);
    virtual double s_dpdrho_i_const_T(const Gas_data &Q, int isp, int &status);
    virtual double s_dpdT_i_const_rho(const Gas_data &Q, int itm, int &status);
    virtual double s_dedT_const_v(const Gas_data &Q, int &status);
    virtual double s_dhdT_const_p(const Gas_data &Q, int &status);
    virtual double s_gas_constant(const Gas_data &Q, int &status);
    virtual double s_mixture_molecular_weight(const Gas_data &Q);
    virtual double s_molecular_weight(int isp);
    virtual double s_internal_energy(const Gas_data &Q, int isp) = 0;
    virtual double s_enthalpy(const Gas_data &Q, int isp) = 0;
    virtual double s_entropy(const Gas_data &Q, int isp) = 0;
    virtual double s_modal_enthalpy(const Gas_data &Q, int isp, int itm);
    virtual double s_modal_Cv(Gas_data &Q, int itm);

private:
    // Local classes.
    class dTdp_functor : public Univariate_functor {
    public:
	dTdp_functor(Gas_model &g, Gas_data &Q)
	    : g_(g), Q_(Q) {}

	double operator()(double p)
	{
	    Q_.p = p;
	    g_.eval_thermo_state_rhop(Q_);
	    return Q_.T[0];
	}
	    
	Gas_model &g_;
	Gas_data &Q_;
    };

    class dTdrho_functor : public Univariate_functor {
    public:
	dTdrho_functor(Gas_model &g, Gas_data &Q)
	    : g_(g), Q_(Q) {}

	double operator()(double rho)
	{
	    Q_.rho = rho; 
	    g_.eval_thermo_state_rhop(Q_);
	    return Q_.T[0];
	}

	Gas_model &g_;
	Gas_data &Q_;
    };

    class dpdrho_functor : public Univariate_functor {
    public:
	dpdrho_functor(Gas_model &g, Gas_data &Q)
	    : g_(g), Q_(Q) {}

	double operator()(double rho)
	{
	    Q_.rho = rho; 
	    g_.eval_thermo_state_rhoT(Q_);
	    return Q_.p;
	}

	Gas_model &g_;
	Gas_data &Q_;
    };

    class dpdrho_i_functor : public Univariate_functor {
    public:
	dpdrho_i_functor(Gas_model &g, Gas_data &Q, int isp)
	    : g_(g), Q_(Q), isp_(isp) {}

	double operator()(double rho_i)
	{
	    Q_.massf[isp_] = rho_i/Q_.rho;
	    g_.eval_thermo_state_rhoT(Q_);
	    return Q_.p;
	}

	Gas_model &g_;
	Gas_data &Q_;
	int isp_;

    };


    class dedT_functor : public Univariate_functor {
    public:
	dedT_functor(Gas_model &g, Gas_data &Q)
	    : g_(g), Q_(Q) {}

	double operator()(double T)
	{
	    Q_.T[0] = T;
	    g_.eval_thermo_state_rhoT(Q_);
	    return Q_.e[0];
	}
	
	Gas_model &g_;
	Gas_data &Q_;
    };

    class dhdT_functor : public Univariate_functor {
    public:
	dhdT_functor(Gas_model &g, Gas_data &Q)
	    : g_(g), Q_(Q) {}
	
	double operator()(double T)
	{
	    Q_.T[0] = T;
	    g_.eval_thermo_state_pT(Q_);
	    return Q_.e[0] + Q_.p/Q_.rho;
	}

	Gas_model &g_;
	Gas_data &Q_;
    };

};

Gas_model* create_gas_model(std::string cfile);
// For swig as delete_Gas_model doesn't work
void call_gas_model_deconstructor( Gas_model *gm );

//int declare_model(lua_State *L);
int declare_species(lua_State *L);
lua_State* initialise_lua_State();
void convert_massf2molef(const std::vector<double> &massf,
			 const std::vector<double> &M,
			 std::vector<double> &molef);
void convert_molef2massf(const std::vector<double> &molef,
			 const std::vector<double> &M, 
			 std::vector<double> &massf);
void convert_massf2conc(double rho,
			const std::vector<double> &massf,
			const std::vector<double> &M,
			std::vector<double> &c);
// This overloaded version is only here because of the
// ODE library.  When the ODE library is changed to use
// vector<double>, then this can be removed.
void convert_massf2conc(double rho,
			const std::vector<double> &massf,
			const std::vector<double> &M,
			std::valarray<double> &c);
void 
convert_conc2massf(double rho,
		   const std::vector<double> &c,
		   const std::vector<double> &M,
		   std::vector<double> &massf);
// Same comment as above.
void convert_conc2massf(double rho,
			const std::valarray<double> &c,
			const std::vector<double> &M,
			std::vector<double> &massf);
void convert_conc2molef(double rho_bar,
			const std::valarray<double> &c,
			std::vector<double> &molef);

// python friendly versions
std::vector<double> convert_massf2molef(const std::vector<double> &massf,
					const std::vector<double> &M);
std::vector<double> convert_molef2massf(const std::vector<double> &molef,
					const std::vector<double> &M);
std::vector<double> convert_massf2conc(double rho,
				       const std::vector<double> &massf,
				       const std::vector<double> &M);
std::vector<double> convert_conc2massf(double rho,
	                               const std::vector<double> &c,
                                       const std::vector<double> &M);

double calculate_molecular_weight(const std::vector<double> &massf,
				  const std::vector<double> &M);
double calculate_gas_constant(const std::vector<double> &massf,
			      const std::vector<double> &M);

#endif
