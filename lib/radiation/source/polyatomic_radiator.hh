/** \file diatomic_radiator.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 20-Feb-2013: Initial version
 *
 *  \brief Declarations for polyatomic radiator classes
 *
 **/

#ifndef POLYATOMIC_RADIATOR_HH
#define POLYATOMIC_RADIATOR_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"
#include "radiator.hh"

class PolyatomicElecLev : public ElecLev {
public:
    /// \brief Constructor
    PolyatomicElecLev( int ilev, double E, int g, lua_State * L );
    
    /// \brief Deconstructor
    virtual ~PolyatomicElecLev() {};

public:
    /// \brief String representation
    std::string string();

    /// \brief Calculate equilibrium total partion function (QelQvQr)
    double calculate_equilibrium_Q_total( double T );

    /// \brief Calculate rovibrational partion function (QvQr)
    double calculate_Q_vib( double T_vib, double T_rot );

    /// \brief Calculate and store rovibrational partion function (QvQr)
    double calculate_and_store_Q_vib( double T_vib, double T_rot );

private:
    /* Energy terms */
    double D;           // Dissociation energy

    /* Partition functions */
    double QvQr;
};

class PolyatomicRadiator : public Radiator {
public:
    /// \brief Constructor
    PolyatomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    virtual ~PolyatomicRadiator();
    
public:
    /// \brief Set the electron index variable
    void set_e_index( int iel );
    
    /// \brief Return pointer to electronic level as an 'ElecLev'
    ElecLev * get_elev_pointer( int ie );

    /// \brief External access for the dissociation energy
    double get_D();

protected:
    /* Energy terms */
    double D;           // Dissociation energy
    
protected:
    /*** Thermodynamics functions ***/
    
    /// \brief Calculate and store internal partition functions (QvQr[], Q_el)
    void calculate_Q_int( Gas_data &Q );
    
    /// \brief Calculate the total equilibrium partition function for this radiator
    double calculate_total_equil_partition_function( double T );
    
    /*** Initiatisation functions ***/
    
    /// \brief Initialise radiation mechanisms
    void initialise_mechanisms( Gas_data &Q );
    
    /*** Emission coefficient functions ***/
    
    /// \brief Calculate the spectrally unresolved emission coeffient
    virtual double calculate_unresolved_emission_coefficient( Gas_data &Q );  
    
    /// \brief Calculate the spectrally unresolved optically variable emission coeffient
    virtual double calculate_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u);

    /*** Spectral generation functions ***/
    
    /// \brief Fill the frequency vector with optimized spectral points
    void spectral_distribution( std::vector<double> &nus );
    
    /// \brief Calculate and store the spectrum for this radiator
    void calculate_spectrum( Gas_data &Q, CoeffSpectra &X ); 
    
    /*** Analysis functions ***/
    
    /// \brief Return string of line-widths for the given gas-state
    std::string line_width_string( Gas_data &Q );
    
public:
    int iTv, iTr;
    
protected:
    // std::vector<PolyatomicSystem*> systems;
    std::vector<PolyatomicElecLev*> elevs;
};

class BoltzLinearPolyatomicRadiator : public PolyatomicRadiator {
public:
    /// \brief Constructor
    BoltzLinearPolyatomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    ~BoltzLinearPolyatomicRadiator();
    
private:
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };
};

#endif
