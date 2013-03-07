// Author: Daniel F. Potter
// Version: 24-Mar-2010
//          Ported from lib/radiation/source/diatomic_radiator.hh

#ifndef DIATOM_ELEC_LEV_HH
#define DIATOM_ELEC_LEV_HH

#define E_VIB_RELATIVE_TO_GROUND_STATE 1

#include <vector>

#include "gas_data.hh"

class Diatom_electronic_level {
public:
    /// \brief Constructor
    Diatom_electronic_level( std::vector<double> lev_data );
    
    /// \brief Deconstructor
    virtual ~Diatom_electronic_level() {};
    
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
    
    /// \brief Calculate anharmonic vibrational energy for given quantum state
    double eval_E_vib( int iV );
    
    /// \brief Calculate rotational energy assuming a singlet state
    double eval_E_rot( int iV, int iJ );
    
    /// \brief Calculate Bv term for the given vibrational level
    double eval_B_v( int iV );
    
    /// \brief Calculate the rotational partition function for this electronic level
    double eval_Qr( double T_rot );

    /// \brief Calculate the rovibrational partition function for this electronic level
    double eval_QvQr( double T_vib, double T_rot );

    /// \brief Calculate the rovibronic partition function for this electronic level
    double eval_QeQvQr( double T_el, double T_vib, double T_rot );

public:
    /* Fundamental level data */
    double N;
    double E;
    int g;
    
    /* Quantum state */
    int spin;		// State spin multiplicity
    int Lambda;		// ND orbital anglular momentum
    double A_spin;	// Spin-angular momentum coupling constant A ('spin-orb')
    
    /* Energy terms */
    double D;
    double r_e;
    
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
};

#endif
