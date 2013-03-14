/** \file atomic_line.hh
 *  \ingroup radiation
 *
 *  \author Daniel F. Potter
 *  \version 21-Aug-07
 #  \brief Declarations for AtomicLine class
 *
 **/

#ifndef ATOMIC_LINE_HH
#define ATOMIC_LINE_HH

#include <string>
#include <vector>

#include "../../nm/source/fobject.hh"

#include "spectra_pieces.hh"

#define ATOMIC_VOIGT_PROFILE_METHOD 1      /* [0] Accurate Whiting expression, [1] approx. Whiting expression */
#define ATOMIC_STARK_WIDTH          0      /* [0] Johnston 2006 curve fit, [1] Cowley 1971, [2] Arnold 1979 */
#define LIMITED_ATOMIC_LINE_EXTENT  1      /* Line extents are unlimited [0] or limited [1]       */

// Forward declaration of AtomicElecLev
class AtomicElecLev;

class AtomicLine {
public:
    /// \brief Default constructor
    AtomicLine() {};
    
    /// \brief Constructor
    AtomicLine( std::vector<double> line_data, double m_w, double I, int npoints, int nwidths, double beta );
    
    /// \brief Deconstructor
    ~AtomicLine();
    
public:
    /// \brief String representation
    std::string string();
    
    /// \brief Initialisation function
    void initialise( double T, double T_e, double p, double N_hvy, double N_e, double mw_av );
    
    /// \brief Fill the frequency vector with optimized spectral points
    void spectral_distribution( std::vector<double> &nus );

    /// \brief Calculate and store the spectra for the given spectral distribution
    void calculate_spectrum( CoeffSpectra &X );
    
    /// \brief Lorentzian line-shape function
    double get_lorentz_point( double delta_nu );
    
    /// \brief Doppler line-shape function
    double get_doppler_point( double delta_nu );
    
    /// \brief Voigt line-shape function
    double get_voigt_point( double delta_nu );
    
    /// \brief List of line-width values for given gas-state
    std::string line_width_string( double T, double T_e, double p, double N_hvy, double N_e, double mw_av );
    
public:
    /// \brief Calculate Lorentzian line-width
    double calculate_lorentz_width( double T, double T_e, double p, double N_hvy, 
    	                            double N_e, double N_l, double mw_av );
    
    /// \brief Calculate resonance line-width 
    double calculate_resonance_width( double N_l );
    
    /// \brief Calculate [Johnannes Diderik] van der Waals line-width 
    double calculate_vanderwaals_width( double T, double p, double mw_av, double N_hvy );
    
    /// \brief Calculate Stark line-width 
    double calculate_Stark_width( double T_e, double N_e );
    
    /// \brief Calculate natural line-width 
    double calculate_natural_width();
    
    /// \brief Calculate Doppler line-width
    double calculate_doppler_width( double T );
    
    /// \brief Calculate Voigt line-width (requires gamma_L and gamma_D)
    double calculate_voigt_width();
    
public:
    /* Minimum data for calculating emissive energy of a line */
    double nu_ul;                        /**< \brief Transition frequency                 */
    double E_l;                          /**< \brief Lower electronic state energy (J)    */
    double E_u;                          /**< \brief Upper electronic state energy (J)    */
    int g_l;                             /**< \brief Lower electronic state degeneracy    */
    int g_u;                             /**< \brief Upper electronic state degeneracy    */
    double A_ul;                         /**< \brief Spon. emission einstein coefficient  */
    double B_lu;                         /**< \brief Induced absorp. einstein coefficient */
    double f_lu;                         /**< \brief Induced absorp. oscillator strength  */
    
    /* Species data for calculating line widths */
    double m_w;
    double I;
    
    /* Data for CR modelling (also for j_ul and kappa_lu calculation) */
    int ie_l, ie_u;
    AtomicElecLev * elev_l;
    AtomicElecLev * elev_u;
    int type;
    
    /* Further data for describing the line-shape */
    double gamma_L;                      /**< \brief Combined Lorentz (half) half-width   */
    double gamma_D;                      /**< \brief Doppler (half) half-width            */
    double gamma_V;                      /**< \brief Voigt (half) half-width              */
    
    /* Line intensity (both emission and absorption coefficients). */
    double j_ul;
    double kappa_lu;

    /* Adaptive spectral grid parameters */
    int npoints;
    int nwidths;
    RobertsClusterFunction * rcf;
};

#endif

