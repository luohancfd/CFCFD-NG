/** \file planck_radiator.hh
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: New radiation2 version
 *
 *  \brief Declarations for planck radiator classes
 *
 **/

#ifndef PLANCK_RADIATOR_HH
#define PLANCK_RADIATOR_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"
#include "radiator.hh"

class PlanckRadiator : public Radiator {
public:
    /// \brief Constructor
    PlanckRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    virtual ~PlanckRadiator();
    
public:
    /// \brief Set the electron index variable
    void set_e_index( int iel );
    
    /// \brief Return pointer to electronic level as an 'ElecLev'
    ElecLev * get_elev_pointer( int ie );
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };
    
private:
    /*** Thermodynamics functions ***/
    
    /// \brief Calculate and store internal partion functions
    void calculate_Q_int( Gas_data &Q );
    
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
    
    /// \brief Calculate the total equilibrium partition function for this radiator
    double calculate_total_equil_partition_function( double T );
    
    /*** Initiatisation functions ***/
    
    /// \brief Initialise radiation mechanisms
    void initialise_mechanisms( Gas_data &Q );
    
    /*** Emission coefficient functions ***/
    
    /// \brief Calculate the spectrally unresolved emission coeffient
    virtual double calculate_unresolved_emission_coefficient( Gas_data &Q ); 
    
    double calculate_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u)
    { return 0.0; }
    
    /*** Spectral generation functions ***/
    
    /// \brief Fill the frequency vector with optimized spectral points
    void spectral_distribution( std::vector<double> &nus );
    
    /// \brief Calculate and store the spectrum
    void calculate_spectrum( Gas_data &Q, CoeffSpectra &X ); 
    
    /*** Analysis functions ***/
    
    /// \brief Return string of line-widths for the given gas-state
    std::string line_width_string( Gas_data &Q );
    
private:
    double kappa_const;
};

#endif
