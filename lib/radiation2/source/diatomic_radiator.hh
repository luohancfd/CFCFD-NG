/** \file diatomic_radiator.hh
 *  \ingroup radiation2
 *
 *  \author Daniel F. Potter
 *  \version 10-Aug-2009: New radiation2 version
 *
 *  \brief Declarations for diatomic radiator classes
 *
 **/

#ifndef DIATOMIC_RADIATOR_HH
#define DIATOMIC_RADIATOR_HH

#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"
#include "radiator.hh"
#include "diatomic_system.hh"
#include "cr_reactions.hh"

/*** Control flags ***/
#define EXACT_Q_ROT 0			// Exact (1) or approx (0) rotational partition functions (NOTE: requires LIMIT_V_TO_GLOBAL_MAX==1)
#define LIMIT_V_TO_GLOBAL_MAX 0		// Limit vibrational quantum numbers to that of the tabulated J_max
#define USE_APPROX_DISSOCIATION_LIMIT 0	// (0) use given dissociation energy, (1) use D ~ omega_e**2 / ( 4*xomega_e ) [(1) is more conservative]
#define MODEL_HUND_AB_TRIPLETS 0

/*** Transition types ***/
#define HUND_A            0
#define SIGMA_TRIPLET     1
#define HUND_B_DOUBLET    2
#define HUND_B_TRIPLET    3
#define HUND_AB_DOUBLET   4
#define HUND_AB_TRIPLET   5

/* For the QSS calculation */
#define ENFORCE_DIATOMIC_BOLTZ_QSS_LIMIT 1

// Minimum j_band / j_system for band inclusion
#define F_BAND_LIMIT 1.0e-3

class DiatomicElecLev : public ElecLev {
public:
    /// \brief Constructor
    DiatomicElecLev( int ilev, std::vector<double> lev_data, bool homonuclear, double I_spin );
    
    /// \brief Deconstructor
    virtual ~DiatomicElecLev() {};
    
public:
    /// \brief Calculate maximum vibrational quantum number
    void calculate_V_max();
    
    /// \brief Test the maximum vibrational quantum number
    bool test_V_max();
    
    /// \brief Evaluate the Morse plus centrifugal potential
    double eval_potential_curve( int iJ, double r );
    
    /// \brief Evaluate the 1st derivative of the Morse plus centrifugal potential wrt r
    double eval_potential_curve_first_derivative( int iJ, double r );
    
    /// \brief Evaluate the 2nd derivative of the Morse plus centrifugal potential wrt r
    double eval_potential_curve_second_derivative( int iJ, double r );
    
    /// \brief Solve for the r_max value of the potential curve
    double solve_for_potential_curve_r_max( int iJ );

    /// \brief Calculate maximum rotational quantum numbers
    void calculate_J_max();
    
    /// \brief Test the maximum rotational quantum number
    bool test_J_max( int iV );
    
    /// \brief Write some potential curves to file
    void write_potential_curves( std::string species_name, int ilev );
    
    /// \brief String representation
    std::string string();
    
    /// \brief Calculate equilibrium total partion function (QelQvQr)
    double calculate_equilibrium_Q_total( double T );
    
    /// \brief Calculate rovibrational partion function (QvQr)
    double calculate_Q_vib( double T_vib, double T_rot );
    
    /// \brief Calculate and store rovibrational partion function (QvQr)
    double calculate_and_store_Q_vib( double T_vib, double T_rot );
    
    /// \brief Calculate anharmonic vibrational energy for given quantum state
    double calculate_E_vib( int iV );
    
    /// \brief Calculate rotational partion function
    double calculate_Q_rot( double T, int iV );
    
    /// \brief Calculate rotational energy assuming a singlet state
    double calculate_E_rot_HundA( int iV, int itJ );
    
    /// \brief Calculate effective rotational energy assuming a triplet state
    double calculate_E_rot_SigmaTriplet( int iV, int itJ );
    
    /// \brief Calculate rotational energy assuming Hund case b
    double calculate_E_rot_HundBDoublet( int iV, int itJ, int itSig );
    
    /// \brief Calculate rotational energy assuming Hund's intermediate a/b case
    double calculate_E_rot_HundABDoublet( int iV, int itJ, int itSig );
    
    /// \brief Calculate rotational energy assuming Hund's intermediate a/b case
    double calculate_E_rot_HundABTriplet( int iV, int itJ, int itSig );
    
    /// \brief Return pre-calculated rovibrational partition function
    double get_QvQr() { return QvQr; }
    
    /// \brief Calculate the Boltzmann population of the given vibrational level
    double calculate_N_vib( double T_vib, double T_rot, int iV );
    
    /// \brief Calculate the rotational Boltzmann population assuming a singlet state
    double calculate_N_rot_HundA( double T_vib, double T_rot, int iV, int itJ, bool apply_LAF=false );
    
    /// \brief Calculate the effective rotational Boltzmann population assuming an (lumped) triplet sigma state
    double calculate_N_rot_SigmaTriplet( double T_vib, double T_rot, int iV, int itJ, bool apply_LAF=false );
    
    /// \brief Calculate the rotational Boltzmann population assuming Hund case b
    double calculate_N_rot_HundBDoublet( double T_vib, double T_rot, int iV, int itJ, int itSig, bool apply_LAF=false );
    
    /// \brief Calculate the rotational Boltzmann population assuming Hund's intermediate a/b case
    double calculate_N_rot_HundABDoublet( double T_vib, double T_rot, int iV, int itJ, int itSig, bool apply_LAF=false );
    
    /// \brief Calculate the rotational Boltzmann population assuming Hund's intermediate a/b case
    double calculate_N_rot_HundABTriplet( double T_vib, double T_rot, int iV, int itJ, int itSig, bool apply_LAF=false );
    
    /// \brief Evaluate line alternation factor
    double line_alternation_factor( int itJ, int itSig=0 );
    
    /// \brief Calculate Bv term for the given vibrational level
    double calculate_B_v( int iV );
    
    /// \brief External access to J_max vector
    int get_J_max( int iV );
    
    /// \brief External access to lambda
    int get_lambda();
    
    /// \brief External access to spin
    int get_spin();
    
    /// \brief External access to A_spin
    double get_A_spin();
    
    /// \brief External access to presence of gamma_V terms
    bool has_gamma_Vs();
    
    /// \brief External access to gamma_V vector
    double get_gamma_V( int iV );
    
    /// \brief External access to D_e
    double get_D_e();
    
    /// \brief External access to V_max
    int get_V_max();
    
    /// \brief Kronecker delta function for lambda splitting
    int eval_kronecker_delta();
    
    /// \brief Homonuclear flag
    bool get_homonuclear();
    
    /// \brief External access to dissociation energy
    double get_D();
    
    /// \brief External access to rotational coupling constant
    double get_B_e();
    
private:
    /* Nuclear symmetry */
    bool homonuclear;
    double I_spin;
    int P_pm;
    int P_gu;
    
    /* Quantum state */
    int spin;		// State spin multiplicity
    int lambda;		// ND orbital anglular momentum
    double A_spin;	// Spin-angular momentum coupling constant A ('spin-orb')
    
    /* Energy terms */
    double D;		// Dissociation energy
    double r_e;		// potential well minimum
    
    /* Klein-Dunham coefficients */
    /* a. Vibrational */
    int V_max;
    double omega_e;
    double xomega_e;
    double yomega_e;
    double zomega_e;
    
    /* b. Rotational */
    std::vector<int> J_max;
    double B_e;
    double alpha_e;
    double D_e;
    double beta_e;
    
    /* c. Spin-rotation constants (for multiplet Sigma states only) */
    std::vector<double> gamma_V;
    
    /* Partition functions */
    double QvQr;
};

class DiatomicRadiator : public Radiator {
public:
    /// \brief Constructor
    DiatomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    virtual ~DiatomicRadiator();
    
public:
    /// \brief Set the electron index variable
    void set_e_index( int iel );
    
    /// \brief Return pointer to electronic level as an 'ElecLev'
    ElecLev * get_elev_pointer( int ie );
    
    /// \brief Return pointer to electronic level as an 'DiatomicElecLev'
    DiatomicElecLev * get_diatomic_elev_pointer( int ie );
    
    /// \brief Return name of a system
    std::string get_system_name( int isys );
    
    /// \brief Return name of a system
    DiatomicSystem * get_system_pointer( int isys );
    
    /// \brief External access for the dissociation energy
    double get_D();
    
    /// \brief Calculate the average wavenumber in cm-1 for the requested vibronic transition
    double calculate_vibronic_wavenumber( int ie_u, int ie_l, int Vu, int Vl );
    
    /// \brief Calculate the average wavelength in nm for the requested system band
    double calculate_vibronic_wavelength( std::string sys_name, int Vu, int Vl );
    
    /// \brief Calculate the Kronecker delta for the requested electronic level
    int calculate_kronecker_delta_for_elev( int ie );
    
    /// \brief Calculate the effective transitive probability for a system
    double calculate_system_transition_probability( int isys, double Tv );
    
protected:
    /*** Initialisation functions ***/
    
    /// \brief Read electronic level data
    void read_elevel_data( lua_State * L );
    
    /// \brief Read radiative system data
    void read_system_data( lua_State * L );
    
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
    
    /// \brief Calculate the spectrally unresolved emission coeffient
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
    double D;
    double I_spin;	// Nuclear spin
    bool homonuclear;
    std::vector<DiatomicSystem*> systems;
    std::vector<DiatomicElecLev*> elevs;
};

class BoltzDiatomicRadiator : public DiatomicRadiator {
public:
    /// \brief Constructor
    BoltzDiatomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    ~BoltzDiatomicRadiator();
    
private:
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string() { return ""; };
};

class QSSDiatomicRadiator : public DiatomicRadiator {
public:
    /// \brief Constructor
    QSSDiatomicRadiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    ~QSSDiatomicRadiator();
    
public:
    /// \brief Set the radiator pointers required by EII mechanisms
    void set_radiator_pointers( std::vector<Radiator*> radiators );
    
    /// \brief Write currently stored level populations to file
    void level_population_file( Gas_data &Q, int index );
    
    /// \brief Calculate Boltzmann population for a level
    double eval_Boltzmann_population_for_level( Gas_data &Q, int ilev );
    
    /// \brief Calculate dissociation equilibrium population for a level
    double eval_DE_population_for_level( Gas_data &Q, int ilev );
    
    /// \brief Collisional-radiative model in LaTeX format
    std::string CR_model_latex_string();
    
private:
    /// \brief Create the CR reactions
    int create_reactions( lua_State * L );
    
    /// \brief Calculate electronic state number densities
    void calculate_n_e( Gas_data &Q );
    
public:
    std::vector<NoneqElecLev*> noneq_elevs;    
    
private:
    Radiator * elec;
    Radiator * atom_A;
    Radiator * atom_B;
    std::vector<CR_Reaction*> reactions;
    
    Valmatrix * dGdy;
    std::valarray<double> C;
    std::valarray<double> y_out;
    
    double T_lower;
};

// Some helper functions
int get_diatomic_transition_type( DiatomicElecLev * elev_u, DiatomicElecLev * elev_l );

#endif
