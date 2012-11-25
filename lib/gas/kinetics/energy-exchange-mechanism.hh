// Author: Daniel F. Potter
// Date: 18-Nov-2009

#ifndef ENERGY_EXCHANGE_MECHANISMS_HH
#define ENERGY_EXCHANGE_MECHANISMS_HH

#include <string>
#include <valarray>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../models/gas_data.hh"
#include "../models/chemical-species.hh"
#include "../models/chemical-species-library.hh"

#include "energy-exchange-relaxation-time.hh"

class Energy_exchange_mechanism {
public:
    /// \brief Default constructor
    Energy_exchange_mechanism();

    /// \brief Default destructor
    virtual ~Energy_exchange_mechanism();

    void compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { tau_ = specific_compute_relaxation_time(Q,molef); }
    
    double compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);

    double py_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);

    double get_tau()
    { return tau_; }

protected:
    double tau_;
    virtual double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef) = 0;
    virtual double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef) = 0;
};

class VT_exchange : public Energy_exchange_mechanism {
public:
    VT_exchange(lua_State *L, int ip, int imode);

    ~VT_exchange();
private:
    int ip_; 	// vibrating species index
    int iT_;	// translational energy mode index
    int iTv_;	// vibrational energy mode index
    Species_energy_mode *p_vib_;	// pointer to species p vibrational energy mode
    Relaxation_time *tau_VT_;
    double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return tau_VT_->compute_relaxation_time(Q, molef); }

    double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
};

// class Polymodal_VT_exchange : public Energy_exchange_mechanism {
// public:
//     Polymodal_VT_exchange( lua_State *L );

//     ~Polymodal_VT_exchange();
// private:
//     int ip_; 	// vibrating species index
//     int iT_;	// translational energy mode index
//     int iTv_;	// vibrational energy mode index
//     std::vector<Species_energy_mode*> p_vib_;	// pointer to species p vibrational energy mode
//     Relaxation_time * tau_VT_;
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_VT_->compute_relaxation_time(Q,molef); }

//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

// class ET_exchange : public Energy_exchange_mechanism {
// public:
//     ET_exchange( lua_State *L );

//     ~ET_exchange();
// private:
//     int iT_;	// translational temperature index
//     int iTe_;	// electronic temperature index
//     int ie_;	// electron index
//     Relaxation_time * tau_ET_;
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_ET_->compute_relaxation_time(Q,molef); }

//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

class VV_THO_exchange : public Energy_exchange_mechanism {
public:
    VV_THO_exchange(lua_State *L, int ip, int imode);
    ~VV_THO_exchange();
private:
    int ip_, iq_, iT_, iTvp_, iTvq_;
    Truncated_harmonic_vibration *p_vib_, *q_vib_;
    double theta_v_p_, theta_v_q_;
    Relaxation_time *tau_VV_;
    
    double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return tau_VV_->compute_relaxation_time(Q,molef); }

    double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
};

// class VV_HO_exchange : public Energy_exchange_mechanism {
// public:
//     VV_HO_exchange( lua_State *L );

//     ~VV_HO_exchange();
// private:
//     // Global data
//     int iT_;					// translational temperature
    
//     // Species p data
//     int ip_;					// molecule p species index
//     int iTvp_;					// molecule p vibratonal temperature index
//     Harmonic_vibration * p_vib_;		// pointer to p vibrational mode
    
//     // Species q data
//     int iTvq_;					// molecule q vibrational temperature index
//     Harmonic_vibration * q_vib_; 		// pointer to p vibrational mode
    
//     Relaxation_time * tau_VV_;
    
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_VV_->compute_relaxation_time(Q,molef); }
    
//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

// class VE_exchange : public Energy_exchange_mechanism {
// public:
//     VE_exchange( lua_State *L );
    
//     ~VE_exchange();
// private:
//     // Electron data
//     int iTe_;
//     int ie_;
    
//     // Vibrator data
//     int iTv_;
//     int iv_;
//     Species_energy_mode * v_vib_;
    
//     Relaxation_time * tau_VE_;
    
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_VE_->compute_relaxation_time(Q,molef); }
    
//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

// class EV_exchange : public Energy_exchange_mechanism {
// public:
//     EV_exchange( lua_State *L );
    
//     ~EV_exchange();
// private:
//     // Electron data
//     int iTe_;
//     int ie_;
    
//     // Vibrator data
//     int iTv_;
//     int iv_;
//     Species_energy_mode * v_vib_;
    
//     Relaxation_time * tau_EV_;
    
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_EV_->compute_relaxation_time(Q,molef); }
    
//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

// class RT_exchange : public Energy_exchange_mechanism {
// public:
//     RT_exchange( lua_State *L );

//     ~RT_exchange();
// private:
//     int ip_; 		// rotating species index
//     int iT_;		// translational energy mode index
//     int iTr_;		// rotational energy mode index
//     double Z_R_inf_;	// reference collision number
//     double theta_r_;	// characteristic temperature for rotation
    
//     Species_energy_mode * p_rot_;	// pointer to species p rotational energy mode
//     Relaxation_time * tau_RT_;
    
//     double calculate_collision_number( double T );
    
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_RT_->compute_relaxation_time(Q,molef); }

//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

// class RE_exchange : public Energy_exchange_mechanism {
// public:
//     RE_exchange( lua_State *L );

//     ~RE_exchange();
// private:
//     int iTr_;	// rotational temperature index
//     int iTe_;	// electronic temperature index
//     int ie_;	// electron index
//     Relaxation_time * tau_RE_;
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_RE_->compute_relaxation_time(Q,molef); }

//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

// class ER_exchange : public Energy_exchange_mechanism {
// public:
//     ER_exchange( lua_State *L );

//     ~ER_exchange();
// private:
//     int iTr_;	// rotational temperature index
//     int iTe_;	// electronic temperature index
//     int ie_;	// electron index
//     Relaxation_time * tau_RE_;
//     double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
//     { return tau_RE_->compute_relaxation_time(Q,molef); }

//     double specific_compute_rate(const std::valarray<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

Energy_exchange_mechanism* create_energy_exhange_mechanism(lua_State *L, int ip, int imode);

// Energy_exchange_mechanism* create_energy_exhange_mechanism_from_file( std::string input_file );

#endif
