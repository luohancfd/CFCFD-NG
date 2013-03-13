/** \file radiator.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 21-Aug-07: Initial implementation
 *           06-Jul-09: Improved port from lib/radiation
 *
 *  \brief Declarations for class describing a generic radiating atom/molecule/electron
 *
 **/

#ifndef RADIATOR_HH
#define RADIATOR_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "../../gas/models/gas_data.hh"

#include "spectra_pieces.hh"
#include "photoionisation.hh"

class ElecLev {
public:
    /// \brief Constructor
    ElecLev( int i_val, double E_val, int g_val );
    
    /// \brief Deconstructor
    virtual ~ElecLev();
    
public:
    /// \brief E access (get) function
    double get_E() { return E; };
    
    /// \brief Q_el access (get) function
    double get_Q_el() { return Q_el; };
    
    /// \brief Q_total access (get) function (needed for diatomics)
    double get_Q_int() { return Q_int; }
    
    /// \brief Q_total access (set) function (needed for diatomics)
    void set_Q_int( double Q_internal ) { Q_int = Q_internal; }
    
    /// \brief N access (set) function
    void set_N( double N_el ) { N = N_el; };
    
    /// \brief N access (get) function
    double get_N() { return N; };
    
    /// \brief g access (set) function
    void set_g( int g_el ) { g = g_el; };
    
    /// \brief g access (get) function
    int get_g() { return g; };
    
    /// \brief String representation
    virtual std::string string();
    
    /// \brief Calculate equilibrium total partion function
    virtual double calculate_equilibrium_Q_total( double T );
    
    /// \brief Electronic partition function
    double calculate_Q_el( double T );
    
    /// \brief Calculate AND store the electronic partition function
    double calculate_and_store_Q_el( double T );
    
    /// \brief Calculate Boltzmann population relative to a pre-calculated level
    double calculate_boltz_N( double T, int g_i, double E_i );
    
    /// \brief Calculate Photo-ionisation cross-section 
    double calculate_sigma_bf( double nu )
    { return PICS_model->eval( nu ); };
    
public:
    /* Fundamental level data */
    int i;
    double N;
    double E;
    int g;
    
    /* Temporary partition function storage */
    double Q_el;
    double Q_int;
    
    /* Photo-ionisation cross-section model */
    PhotoIonisationCrossSectionModel * PICS_model;
};

class NoneqElecLev {
public:
    /// \brief Constructor
    NoneqElecLev( int ilev, int ne_ilev, std::string label, ElecLev * elev, std::vector<ElecLev*> eq_elevs );
    
    /// \brief Deconstructor
    ~NoneqElecLev();
    
public:
    int ilev;
    int ne_ilev;
    std::string label;
    ElecLev * elev;
    std::vector<ElecLev*> eq_elevs;
};

class Radiator {
public:
    /// \brief Constructor
    Radiator( lua_State * L, std::string name );
    
    /// \brief Deconstructor
    virtual ~Radiator();
    
public:
    /*** External access functions ***/
    /// \brief Set the electron index variable
    virtual void set_e_index( int iel );
    
    /// \brief Calculate the total equilibrium partition function for this radiator
    double calc_total_equil_partition_function( double T )
    { return calculate_total_equil_partition_function(T); }
    
    /// \brief Evaluate the translational partition function for this radiator
    double eval_translational_partition_function(Gas_data &Q);
    
    /// \brief Evaluate the translational partition function for this radiator
    double eval_translational_partition_function_from_T( double T );
    
    /// \brief Calculate and store internal partion functions
    void calc_partition_functions(Gas_data &Q)
    { calculate_Q_int(Q); }
    
    /// \brief Calculate and store electronic state number densities
    void calc_elec_pops(Gas_data &Q)
    { calculate_n_e(Q); }
    
    /// \brief Calculate and store the spectrum for this radiator
    void calc_spectrum( Gas_data &Q, CoeffSpectra &X )
    { calculate_spectrum( Q, X ); }
    
    /// \brief Calculate the spectrally unresolved emission coeffient
    double calc_unresolved_emission_coefficient( Gas_data &Q)
    { return calculate_unresolved_emission_coefficient( Q ); }
    
    /// \brief Calculate the spectrally unresolved emission coeffient
    double calc_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u)
    { return calculate_unresolved_OV_emission_coefficient( Q, wavel_switch, Lambda_l, Lambda_u ); }
    
    /// \brief Initialise radiation mechanisms for this radiator
    void init_mechanisms( Gas_data &Q )
    { initialise_mechanisms( Q ); }
    
    /// \brief Return pointer to electronic level as an 'ElecLev'
    virtual ElecLev * get_elev_pointer( int ie ) = 0;
    
    /// \brief Return string of line-widths for the given gas-state
    std::string get_line_width_string( Gas_data &Q )
    { return line_width_string( Q ); }
    
    /// \brief Fill the frequency vector with optimized spectral points
    void get_spectral_distribution( std::vector<double> &nus )
    { return spectral_distribution( nus ); }
    
    void write_level_populations_to_file( Gas_data &Q, int index=0 )
    { return level_population_file(Q,index); }
    
    void prep_spatial_population_file()
    { return prep_x_pop_file(); }
    
    void append_current_populations( double x )
    { return append_current_pops( x ); }
    
    double get_Q_int() { return Q_int; }
    
    bool test_for_optically_allowed_transition( int ilev_i, int ilev_f )
    { return optically_allowed_transition_test( ilev_i, ilev_f ); }
    
    std::string get_CR_model_latex_string()
    { return CR_model_latex_string(); }
    
    double sum_level_populations();

protected:
    /*** Data-input functions ***/
    
    /// \brief Read in the photoionization data for this radiator
    void read_photoionization_data( lua_State * L );
    
    /*** Thermodynamics functions ***/
    
    /// \brief Calculate and store internal partion functions
    virtual void calculate_Q_int(Gas_data &Q) = 0;

    /// \brief Calculate and store electronic state number densities
    virtual void calculate_n_e(Gas_data &Q) = 0;

    /// \brief Calculate the total equilibrium partition function for this radiator
    virtual double calculate_total_equil_partition_function( double T ) = 0;
    
    /*** Initiatisation functions ***/
    
    /// \brief Initialise radiation mechanisms for this radiator
    virtual void initialise_mechanisms( Gas_data &Q ) = 0;
    
    /*** Emission coefficient functions ***/
    
    /// \brief Calculate the spectrally unresolved emission coeffient
    virtual double calculate_unresolved_emission_coefficient( Gas_data &Q ) = 0;   
    
    /// \brief Calculate the spectrally unresolved optically variable emission coeffient
    virtual double calculate_unresolved_OV_emission_coefficient( Gas_data &Q, double wavel_switch, double Lambda_l, double Lambda_u) = 0;
    
    /*** Spectral generation functions ***/
    
    /// \brief Fill the frequency vector with optimized spectral points
    virtual void spectral_distribution( std::vector<double> &nus ) = 0;
    
    /// \brief Calculate and store the spectrum for this radiator
    virtual void calculate_spectrum( Gas_data &Q, CoeffSpectra &X ) = 0;  
    
    /*** Analysis functions ***/
    
    /// \brief Return string of line-widths for the given gas-state
    virtual std::string line_width_string( Gas_data &Q ) = 0;
    
    /// \brief Write currently stored level populations to file
    virtual void level_population_file( Gas_data &Q, int index );
    
    /// \brief Prepare the spatial population file
    void prep_x_pop_file();
    
    /// \brief Append data to the spatial population file
    void append_current_pops( double x );
    
    /// \brief Test for an optically allowed transition
    virtual bool optically_allowed_transition_test( int ilev_i, int ilev_f );
    
    /// \brief Collisional-radiative model in LaTeX format
    virtual std::string CR_model_latex_string() = 0;
    
public:
    std::string name;
    std::string type;
    std::string EPM;
    int isp;
    double m_w;
    double h_f;
    double I;
    int Z;
    int iT;
    int iTe;
    int nlevs;
    int nsys;
    int e_index;
    
    double Q_int;
    double Q_el;
    
    bool bf_ion_flag;
};

Radiator *create_new_radiator( lua_State * L, const std::string name );

Radiator * get_radiator_pointer_from_name( std::vector<Radiator*> radiators, std::string name );

#endif
