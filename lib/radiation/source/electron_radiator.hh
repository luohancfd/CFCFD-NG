/** \file electron_radiator.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: Improved port from old lib/radiation
 *
 *  \brief Declarations for electron radiator classes
 *
 **/

#ifndef ELECTRON_RADIATOR_HH
#define ELECTRON_RADIATOR_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"
#include "radiator.hh"

class ContinuumMechanism {
public:
    /// \brief Constructor
    ContinuumMechanism();
    
    /// \brief Deconstructor
    virtual ~ContinuumMechanism();
    
public:
    virtual void calculate_spectrum( Gas_data &Q, CoeffSpectra &X, double n_elecs ) = 0;
    
    virtual void initialise_mechanisms( Gas_data &Q ) = 0;
};

class FreeFreeHydrogenic : public ContinuumMechanism {
public:
    /// \brief Constructor
    FreeFreeHydrogenic( Radiator * R_ion );
    
    /// \brief Deconstructor
    ~FreeFreeHydrogenic();
    
public:
    void calculate_spectrum( Gas_data &Q, CoeffSpectra &X, double n_elecs );
    
    void initialise_mechanisms( Gas_data &Q );
    
private:
    Radiator * Ri;
};

class BoundFreeHydrogenic : public ContinuumMechanism {
public:
    /// \brief Constructor
    BoundFreeHydrogenic( Radiator * R_ion, Radiator * R_neutral );
    
    /// \brief Deconstructor
    ~BoundFreeHydrogenic();
    
public:
    void calculate_spectrum( Gas_data &Q, CoeffSpectra &X, double n_elecs );
    
    void initialise_mechanisms( Gas_data &Q );
    
private:
    Radiator * Ri;
    Radiator * Rn;
};

class ElectronRadiator : public Radiator {
public:
    /// \brief Constructor
    ElectronRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    ~ElectronRadiator();
    
public:
    /// \brief Return pointer to electronic level as an 'ElecLev'
    ElecLev * get_elev_pointer( int ie );
    
    /// \brief Create the continuum mechanisms from the initialised radiator set
    void create_continuum_mechanisms( std::vector<Radiator*> &radiators );
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };
    
private:
    /*** Thermodynamics functions ***/
    
    /// \brief Calculate and store internal partion function
    void calculate_Q_int( Gas_data &Q );
    
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
    
    /// \brief Calculate the total equilibrium partition function for this radiator
    double calculate_total_equil_partition_function( double T );
    
    /*** Spectral generation functions ***/
    
    /// \brief Fill the frequency vector with optimized spectral points
    void spectral_distribution( std::vector<double> &nus );
    
    /// \brief Calculate and store the spectrum
    void calculate_spectrum( Gas_data &Q, CoeffSpectra &X );
    
    /*** Emission coefficient functions ***/
    
    /// \brief Calculate the spectrally unresolved emission coeffient
    virtual double calculate_unresolved_emission_coefficient( Gas_data &Q ); 
    
    /// \brief Calculate the spectrally unresolved optically variable emission coeffient
    virtual double calculate_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u)
    { return 0.0; }
    
    /*** Initiatisation functions ***/
    
    /// \brief Initialise radiation mechanisms
    void initialise_mechanisms( Gas_data &Q );
    
    /*** Analysis functions ***/
    
    /// \brief Return string of line-widths for the given gas-state
    std::string line_width_string( Gas_data &Q );
    
private:
    ElecLev * elev;
    
    std::vector<std::string> systems_list;
    
    std::vector<ContinuumMechanism*> continuum_mechanisms;
};

#endif
