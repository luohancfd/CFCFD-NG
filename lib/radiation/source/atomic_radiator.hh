/** \file atomic_radiator.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: Improved port from old lib/radiation
 *
 *  \brief Declarations for atomic radiator classes
 *
 **/

#ifndef ATOMIC_RADIATOR_HH
#define ATOMIC_RADIATOR_HH

#include <string>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"

#include "radiator.hh"
#include "atomic_line.hh"
#include "cr_reactions.hh"

// For creating atomic radiative transitions
#define MIN_TRANSITION_PROBABILITY 0.0
#define NO_DATA -1
#define ALLOWED 0
#define FORBIDDEN 1

// For the QSS calculation
#define ENFORCE_ATOMIC_BOLTZ_QSS_LIMIT 1

class AtomicElecLev : public ElecLev {
public:
    /// \brief Constructor
    AtomicElecLev( int ilev, std::vector<double> lev_data );
    
    /// \brief Copy constructor
    AtomicElecLev( AtomicElecLev * lev );
    
    /// \brief Deconstructor
    virtual ~AtomicElecLev() {};
    
public:
    /// \brief String representation
    std::string string();
    
    int get_n() { return n; }
    int get_l() { return l; }
    int get_L() { return L; }
    int get_S() { return S; }
    int get_parity() { return parity; }
    
private:
    /* Quantum state */
    int n;
    int l;
    int L;
    int S;
    int parity;
};

class AtomicRadiator : public Radiator {
public:
    /// \brief Constructor
    AtomicRadiator( lua_State * L, std::string name );

    /// \brief Deconstructor
    virtual ~AtomicRadiator();
    
public:
    /// \brief Set the electron index variable
    void set_e_index( int iel );
    
    /// \brief Return pointer to electronic level as an 'ElecLev'
    ElecLev * get_elev_pointer( int ie );
    
    /// \brief Search for lines in a particular range
    void find_lines_in_spectral_range( double lambda_min, double lambda_max );
    
protected:
    /*** Initialisation functions ***/
    
    /// \brief Read electronic level data
    void read_elevel_data( lua_State * L );
    
    /// \brief Read line data
    void read_line_data( lua_State * L );
    
    /// \brief Find nearest grouped electronic level for E_i
    int find_grouped_E_level( double E_i );
    
    /*** Thermodynamics functions ***/
    
    /// \brief Calculate and store internal partition function (Q_el)
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
    
    /// \brief Calculate and store the spectrum
    void calculate_spectrum( Gas_data &Q, CoeffSpectra &X ); 
    
    /*** Analysis functions ***/
    
    /// \brief Return string of line-widths for the given gas-state
    std::string line_width_string( Gas_data &Q );
    
    /// \brief Determine if a transition is optically allowed
    bool optically_allowed_transition_test( int ilev_i, int ilev_f );
    
protected:
    int nlines;
    
    std::vector<AtomicElecLev*> elevs;
    std::vector<AtomicLine*> lines;
};

class BoltzAtomicRadiator : public AtomicRadiator {
public:
    /// \brief Constructor
    BoltzAtomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    ~BoltzAtomicRadiator();
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };
    
private:
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
};

class QSSAtomicRadiator : public AtomicRadiator {
public:
    /// \brief Constructor
    QSSAtomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    ~QSSAtomicRadiator();
    
public:
    /// \brief Set the radiator pointers required by EII mechanisms
    void set_radiator_pointers( std::vector<Radiator*> radiators );
    
    /// \brief Write currently stored level populations to file
    void level_population_file( Gas_data &Q, int index );
    
    /// \brief Calculate Boltzmann population for a level
    double eval_Boltzmann_population_for_level( Gas_data &Q, int ilev );
    
    /// \brief Calculate Saha population for a level
    double eval_Saha_population_for_level( Gas_data &Q, int ilev );
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };
    
private:
    /// \brief Create the electron impact excitation reactions
    int create_electron_impact_excitation_reactions( lua_State * L, std::string model );
    
    /// \brief Create the electron impact ionization reactions
    int create_electron_impact_ionization_reactions( lua_State * L, std::string model );
    
    /// \brief Create the radiative transition reactions
    int create_radiative_transition_reactions( lua_State * L, std::string model );
    
private:
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
    
private:
    Radiator * ion;
    Radiator * elec;
    std::vector<NoneqElecLev*> noneq_elevs;
    std::vector<CR_Reaction*> reactions;
    
    Valmatrix * dGdy;
    std::valarray<double> C;
    std::valarray<double> y_out;
    
    double T_lower;
};

class FirstOrderLTNEAtomicRadiator : public AtomicRadiator {
public:
    /// \brief Constructor
    FirstOrderLTNEAtomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    ~FirstOrderLTNEAtomicRadiator();
    
public:
    /// \brief Set the ion and electron pointers
    void set_radiator_pointers( std::vector<Radiator*> radiators );
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };
    
private:
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
    
private:
    int nlevs_boltz;
    Radiator * ion;
    Radiator * elec;
};

class NoneqAtomicRadiator : public AtomicRadiator {
public:
    /// \brief Constructor
    NoneqAtomicRadiator( lua_State * L, std::string name );

    /// \brief Deconstructor
    ~NoneqAtomicRadiator();

    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };

private:
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );

private:
    std::vector<int> isp_list;
};

int get_atomic_transition_type( AtomicElecLev * lev_i, AtomicElecLev * lev_f );

#endif
