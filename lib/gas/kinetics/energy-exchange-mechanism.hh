// Author: Daniel F. Potter
// Date: 18-Nov-2009

#ifndef ENERGY_EXCHANGE_MECHANISMS_HH
#define ENERGY_EXCHANGE_MECHANISMS_HH

#include <string>
#include <vector>

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
    
    double compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);

    double py_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);

    double get_tau()
    { return tau_; }

protected:
    double tau_;
    virtual double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef) = 0;
    virtual double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef) = 0;
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

    double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
};

class Polyatomic_VT_exchange : public Energy_exchange_mechanism {
public:
    Polyatomic_VT_exchange(lua_State *L, int ip, int imode);

    ~Polyatomic_VT_exchange();
private:
    int ip_; 	// vibrating species index
    int iT_;	// translational energy mode index
    int iTv_;	// vibrational energy mode index
    std::vector<Species_energy_mode*> p_vib_;	// pointer to species p vibrational energy mode
    Relaxation_time * tau_VT_;
    double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return tau_VT_->compute_relaxation_time(Q, molef); }

    double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
};

class ET_exchange : public Energy_exchange_mechanism {
public:
    ET_exchange( lua_State *L, int ie, int iTe );

    ~ET_exchange();
private:
    int ie_;	// electron index
    int iTe_;	// electronic temperature index
    int iT_;	// translational temperature index
    Relaxation_time * tau_ET_;
    double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return tau_ET_->compute_relaxation_time(Q,molef); }

    double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
};

class ER_exchange : public Energy_exchange_mechanism {
public:
    ER_exchange( lua_State *L, int ie, int iTe );

    ~ER_exchange();
private:
    int ie_;    // electron index
    int iTe_;   // electronic temperature index
    int iTr_;   // rotational temperature index
    Relaxation_time * tau_ER_;
    double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return tau_ER_->compute_relaxation_time(Q,molef); }

    double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
};

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

    double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
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
    
//     double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
// };

class VE_exchange : public Energy_exchange_mechanism {
public:
    VE_exchange(lua_State *L, int ip, int imode);

    ~VE_exchange();
private:
    int ie_;    // electron species index
    int ip_;    // vibrating species index
    int iTe_;   // electron energy mode index
    int iTv_;   // vibrational energy mode index
    Species_energy_mode *p_vib_;        // pointer to species p vibrational energy mode
    Relaxation_time *tau_VE_;
    double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return tau_VE_->compute_relaxation_time(Q, molef); }

    double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
};

class EV_exchange : public Energy_exchange_mechanism {
public:
    EV_exchange(lua_State *L, int ip, int imode);

    ~EV_exchange();
private:
    int ie_;    // electron species index
    int iq_;    // vibrating species index
    int iTe_;   // electron energy mode index
    int iTv_;   // vibrational energy mode index
    Species_energy_mode *q_vib_;        // pointer to species p vibrational energy mode
    Relaxation_time *tau_EV_;
    double specific_compute_relaxation_time(Gas_data &Q, std::vector<double> &molef)
    { return tau_EV_->compute_relaxation_time(Q, molef); }

    double specific_compute_rate(const std::vector<double> &y, Gas_data &Q, std::vector<double> &molef);
};

Energy_exchange_mechanism* create_energy_exchange_mechanism(lua_State *L, int imode);
Energy_exchange_mechanism* get_mech_from_file(int imech, std::string cfile, Gas_model &g);

#endif
