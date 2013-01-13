/** \file diatomic_system.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 21-Aug-07: Initial version
 *           10-Aug-09: Photaura version
 *  \brief Declarations for DiatomicSystem class
 *
 **/

#ifndef DIATOMIC_SYSTEM_HH
#define DIATOMIC_SYSTEM_HH

#include <string>
#include <vector>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

#include "diatomic_band.hh"

#define TEST_BAND_NU_AV     0
#define WRITE_LINES_TO_FILE 0

// Forward declaration of DiatomicElecLev
class DiatomicElecLev;

class DiatomicSystem {
public:
    /// \brief Constructor
    DiatomicSystem( lua_State * L, std::string name, std::string band_method, 
    		    DiatomicElecLev * elev_l, DiatomicElecLev * elev_u, 
    		    int iT, int iTe, int iTv, int iTr, double m_w, double I  );
    
    /// \brief Deconstructor
    virtual ~DiatomicSystem();
    
public:
    /// \brief Read vibro-electronic transition moment for given vibrational band
    double get_Re_vib_from_file( lua_State *L, int iVu, int iVl, std::string format );
    
    /// \brief Access into 1D band vector
    DiatomicBand * band_pointer( int Vu, int Vl );
    
    /// \brief Initialise the system for spectrum calculation
    void initialize( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    
    /// \brief Calculate the unresolved system emission coefficient
    double calculate_j_ul( double T_vib, double T_rot );
    
    /// \brief Calculate the unresolved optically variable system emission coefficient
    double calculate_OV_j_ul( double T_vib, double T_rot, double wavel_switch, double Lambda_l, double Lambda_u );
    
    /// \brief Calculate the emission and absorption spectra for this system
    void calculate_spectrum( CoeffSpectra &X, double T_vib, double T_rot, double j_av ); 
    
    /// \brief Get the line-width string for this system
    std::string line_width_string( double T, double Te, double p, double N_hvy, double N_elecs, double mw_av );
    
    /// \brief Convert vibrational Einstein coefficient to electronic (vibrational) transition moments
    double convert_VEC_to_ETM( double A, int iVu, int iVl );
    
    /// \brief Convert absorption oscillator strength to electronic (vibrational) transition moments
    double convert_AOS_to_ETM( double q, int iVu, int iVl );
    
    /// \brief Calculate the average frequency of a vibrational band
    double calculate_band_average_frequency( int iVu, int iVl );
    
    /// \brief Calculate the total transition probability for this system
    double calculate_transition_probability( double Tv );
    
public:
    std::string name;
    
    /* Transition data */
    int ie_l, ie_u;
    DiatomicElecLev * elev_u;
    DiatomicElecLev * elev_l;
    int transition_type;
    
    /* Band dimensions */
    int uRe_dim;
    int lRe_dim;
    
    /* Band data */
    std::vector<DiatomicBand*> bands;
    
    /* Other data */
    int iT;
    int iTe;
    int iTv;
    int iTr;
    int e_index;
};

DiatomicSystem * create_new_diatomic_system( lua_State * L, std::string name, 
					     DiatomicElecLev * elev_u, DiatomicElecLev * elev_l,
					     int iT, int iTe, int iTv, int iTr, double m_w, double I );

#endif
