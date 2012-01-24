// Author: Daniel F Potter
// Date: 07-Apr-2010
// Place: Dutton Park, QLD, Oz

#ifndef CHEMISTRY_ENERGY_COUPLING_HH
#define CHEMISTRY_ENERGY_COUPLING_HH

#include <string>
#include <valarray>
#include <vector>

#include "../models/gas_data.hh"

#include "coupling-component.hh"

class Chemistry_energy_coupling {
public:
    /// \brief Normal constructor
    Chemistry_energy_coupling( int isp, std::string mode, std::vector<Coupling_component*> &ccs );
    
    /// \brief Default destructor
    ~Chemistry_energy_coupling();
    
    /// \brief Set the energy and number density from the beginning of the timestep
    void set_e_and_N_old( Gas_data &Q, std::vector<double> &c_old );
    
    /// \brief Compute the new average energy, e_new = e_old + delta_E / N_new
    int update_energy( Gas_data &Q, std::valarray<double> &delta_c, std::vector<double> &c_new );
    
    /// \brief Compute the source term
    double eval_source_term( Gas_data &Q, std::valarray<double> &dcdt );
    
    int get_isp() { return isp_; }
    
    int get_imode() { return imode_; }
    
private:
    int isp_;
    std::string mode_;
    // NOTE: a vector of modes is required to allow for polyatomic vibrational modes
    std::vector<Species_energy_mode*> sems_;
    int imode_;
    double m_;
    double e_old_;
    
    std::vector<Coupling_component*> components_;
};

void create_Chemistry_energy_coupling_for_species_mode( int isp, std::string mode, 
    std::vector<Coupling_component*> &ccs, std::vector<Chemistry_energy_coupling*> &cecs );

#endif
