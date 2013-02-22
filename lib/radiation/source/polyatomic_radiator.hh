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

class PolyatomicVibState {
public:
  PolyatomicVibState( std::string type, int iV, lua_State * L  );

    /// \brief Deconstructor
    virtual ~PolyatomicVibState() {};

public:
    double get_E_vib() { return G_v; }

    int get_J_max() { return J_max; }

    virtual double calculate_E_rot( int iJ, int iK ) = 0;

    double calculate_Q_rot( double T );

protected:
    /* Designation */
    std::string type;
    int iV;
    std::string label;

    /* Spectroscopic data */
    int g;              // vibrational state degeneracy
    double G_v;         // vibrational energy (J)
    double A_v;         // First rotational (inertia) constant (J)
    double B_v;         // Second rotational (inertia) constant (J)
    double C_v;         // Third rotational (inertia) constant (J)
    double D_v;         // First anharmonicity constant (J)
    double H_v;         // Second anharmonicity constant (J)
    int J_max;          // Maximum rotational quantum number
};

class SphericalTopPolyatomicVibState : public PolyatomicVibState {
public:
    SphericalTopPolyatomicVibState( int iV, lua_State * L  );

    /// \brief Deconstructor
    virtual ~SphericalTopPolyatomicVibState() {};

public:
    double calculate_E_rot( int iJ, int iK ); 
};

class SymmetricalTopPolyatomicVibState : public PolyatomicVibState {
public:
    SymmetricalTopPolyatomicVibState( int iV, lua_State * L  );

    /// \brief Deconstructor
    virtual ~SymmetricalTopPolyatomicVibState() {};

public:
    double calculate_E_rot( int iJ, int iK );
};

class AxisymmetricTopPolyatomicVibState : public PolyatomicVibState {
public:
    AxisymmetricTopPolyatomicVibState( int iV, lua_State * L  );

    /// \brief Deconstructor
    virtual ~AxisymmetricTopPolyatomicVibState() {};

public:
    double calculate_E_rot( int iJ, int iK );
};

PolyatomicVibState * create_new_polyatomic_vib_state( std::string type, int iV, lua_State * L );

class PolyatomicElecLev : public ElecLev {
public:
    /// \brief Constructor
    PolyatomicElecLev( int ilev, double E, int g, lua_State * L );
    
    /// \brief Deconstructor
    virtual ~PolyatomicElecLev();

public:
    /// \brief String representation
    virtual std::string string();

    /// \brief Calculate equilibrium total partition function (QelQvQr)
    double calculate_equilibrium_Q_total( double T );

    /// \brief Calculate rovibrational partition function (QvQr)
    double calculate_Q_vib( double T_vib, double T_rot );

    /// \brief Calculate and store rovibrational partition function (QvQr)
    double calculate_and_store_Q_vib( double T_vib, double T_rot );

private:
    /* Symmetry type */
    std::string type;

    /* Energy terms */
    double D;           // Dissociation energy

    /* Vibrational states */
    std::vector<PolyatomicVibState*> vstates;

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
    /*** Initialisation functions ***/

    /// \brief Read electronic level data
    void read_elevel_data( lua_State * L );

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
    /* Energy terms */
    double D;           // Dissociation energy

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
