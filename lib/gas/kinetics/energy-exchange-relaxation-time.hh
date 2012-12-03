// Author: Daniel F. Potter
// Date: 18-Nov-2009

#ifndef ENERGY_EXCHANGE_RELAXATION_TIME_HH
#define ENERGY_EXCHANGE_RELAXATION_TIME_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"

class Relaxation_time {
public:
    Relaxation_time();
    
    virtual ~Relaxation_time();

    double compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return specific_relaxation_time(Q,molef); }

    double compute_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return specific_transition_probability(Q,molef); }
protected:
    virtual double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef) = 0;
    virtual double specific_transition_probability(Gas_data &Q, std::vector<double> &molef) = 0;
};

class VT_MillikanWhite : public Relaxation_time {
public:
    VT_MillikanWhite(lua_State *L);

    ~VT_MillikanWhite();
private:
    double M_p_;
    int ip_;
    double theta_;
    int iT_;
    
    std::vector<bool> homogenous_flags_;
    std::vector<int> iqs_;
    std::vector<double> M_qs_;
    std::vector<double> mus_;
    std::vector<double> a_;
    std::vector<double> b_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class VT_high_temperature_cross_section {
public:
    VT_high_temperature_cross_section() {}

    virtual ~VT_high_temperature_cross_section() {}
    
    double eval_sigma( double T )
    { return s_eval_sigma(T); }
protected:
    virtual double s_eval_sigma( double T ) = 0;
};

class Park_VT_HTCS : public VT_high_temperature_cross_section {
public:
    Park_VT_HTCS(lua_State *L);
    
    ~Park_VT_HTCS();
private:
    double sigma_dash_;
    
    double s_eval_sigma( double T );
};

class Fujita_VT_HTCS : public VT_high_temperature_cross_section {
public:
    Fujita_VT_HTCS();
    
    ~Fujita_VT_HTCS();
private:
    double s_eval_sigma( double T );
};

VT_high_temperature_cross_section *
create_VT_high_temperature_cross_section_model( lua_State *L );

class VT_MillikanWhite_HTC : public Relaxation_time {
public:
    VT_MillikanWhite_HTC(lua_State *L);

    ~VT_MillikanWhite_HTC();
private:
    double M_p_;
    int ip_;
    double theta_;
    int iT_;
    
    std::vector<bool> homogenous_flags_;
    std::vector<int> iqs_;
    std::vector<double> M_qs_;
    std::vector<double> mus_;
    std::vector<double> a_;
    std::vector<double> b_;
    
    VT_high_temperature_cross_section * HTC_model_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ET_AppletonBray : public Relaxation_time {
public:
    ET_AppletonBray(lua_State *L);
    
    ~ET_AppletonBray();
private:
    std::vector<Relaxation_time*> tau_ETs_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ET_Ion : public Relaxation_time {
public:
    ET_Ion(lua_State *L);
    
    ~ET_Ion();
private:
    int ic_;
    int ie_;
    int iTe_;
    double M_c_;
    double M_e_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class ET_Neutral : public Relaxation_time {
public:
    ET_Neutral(lua_State *L);
    
    ~ET_Neutral();
private:
    int ic_;
    double M_c_;
    int iTe_;
    std::vector<double> LT_C_, HT_C_;
    double T_switch_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

// Helper functions for SSH calculations
double collider_distance_A(Diatomic_species &p, Diatomic_species &q);
double collider_distance_B(Diatomic_species &p, Diatomic_species &n);
double xi_correction(Diatomic_species &p, Diatomic_species &n);
double potential_well_A(Diatomic_species &p, Diatomic_species &q);
double potential_well_B(Diatomic_species &p, Diatomic_species &n);
double collision_frequency(double sigma, double mu, double T, double nd);
double SSH_beta(double eps0, double mu, double r0,
		double del_E, double T);
double SSH_D_star(double beta);
double SSH_r_c_star(double beta, double r0);
double SSH_r_c_star2(double D_star, double r0);
double SSH_del_star(double beta, double r0);
double SSH_del_star2(double D_star, double r0);
double SSH_alpha_pq(double mu, double del_E, double del_star);
double SSH_chi_pq(double alpha_pq, double T);
double SSH_Z_0(double del_star, double r_eq);
double SSH_Z_V(double f_m, double mu_pp, double mu_pq, double alpha_pq,
	       double theta_v, double del_E);
double SSH_Z_T(double del_E, double alpha_pq, double T);
double SSH_Z_plus(double eps0, double chi, double T);
double SSH_A_factor(double r_c, double sigma);

class VV_SSH : public Relaxation_time {
public:
    VV_SSH(lua_State *L);
    
    ~VV_SSH();
private:
    // Global data
    int iT_;			// translational temperature
    double sigma_;
    double eps_;
    double mu_;
    double r_;
    double delta_E_;
    
    // Species p data
    int ip_;			// species index
    int iTvp_;			// vibratonal temperature index
    double theta_v_p_;		// characteristic vibrational temperature
    double M_p_;		// molecular weight
    double r_eq_p_;		// equilibrium intermolecular distance (m)
    double f_m_p_;		// Thivet's mass factor
    double mu_p_;		// reduced constituent atom mass
    
    // Species q data
    int iq_;
    int iTvq_;
    double theta_v_q_;
    double M_q_;
    double r_eq_q_;
    double f_m_q_;
    double mu_q_;
    
    double VV_SSH_transition_probability( double T );
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class VV_Candler : public Relaxation_time {
public:
    VV_Candler(lua_State *L);
    
    ~VV_Candler();
private:
    // Global data
    int iT_;			// translational temperature
    double sigma_;		// collision cross-section (dxd)
    double mu_;			// reduced molecular weight
    double P_;			// transition probability
    
    // Species p data
    int ip_;			// species index
    
    // Species q data
    int iq_;			// species index
    double M_q_;		// molecular weight
      
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class VE_Lee : public Relaxation_time {
public:
    VE_Lee(lua_State *L);
    
    ~VE_Lee();
private:
    // Electron data
    int ie_;
    int iTe_;
    
    // Vibrational species data
    int iv_;
    
    std::vector<double> T_switches_;
    std::vector< std::vector<double> > ptau_coeffs_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class RT_Parker : public Relaxation_time {
public:
    RT_Parker(lua_State *L);
    
    ~RT_Parker();
private:
    // Global data
    int iT_;					// translational temperature
    std::vector<double> sigmas_;		// rotational collision cross-sections
    std::vector<double> mus_;			// reduced molecular weights

    // Rotator data
    int ip_;					// species index
    
    // Collider data
    std::vector<int> iqs_;			// species indices
    std::vector<double> M_qs_;			// molecular weights

    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class RE_Abe : public Relaxation_time {
public:
    RE_Abe(lua_State *L);
    
    ~RE_Abe();
private:
    std::vector<Relaxation_time*> tau_REs_;
    
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class RE_Ion : public Relaxation_time {
public:
    RE_Ion(lua_State *L);
    
    ~RE_Ion();
private:
    int ic_;
    int ie_;
    int iTe_;
    double M_c_;
    double M_e_;
    double g_rot_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

class RE_Neutral : public Relaxation_time {
public:
    RE_Neutral(lua_State *L);
    
    ~RE_Neutral();
private:
    int ic_;
    double M_c_;
    int iTe_;
    std::vector<double> C_;
    double g_rot_;
    double specific_relaxation_time(Gas_data &Q, std::vector<double> &molef);
    double specific_transition_probability(Gas_data &Q, std::vector<double> &molef)
    { return 0.0; }
};

Relaxation_time* create_new_relaxation_time( lua_State * L );

#endif
