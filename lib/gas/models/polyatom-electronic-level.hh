// Author: Daniel F. Potter
// Version: 27-July-2012
//          Ported from lib/gas/models/diatom-electronic-level.hh

#ifndef POLYATOM_ELEC_LEV_HH
#define POLYATOM_ELEC_LEV_HH

#include <vector>

#include "gas_data.hh"

class Polyatom_vibrational_mode {
public:
    /// \brief Constructor
    Polyatom_vibrational_mode( int g, double omega, int Vmax );

    /// \brief Deconstructor
    virtual ~Polyatom_vibrational_mode() {};

public:
    /// \brief compute the statistical weight for the vibrational state
    int calculate_p( int iV );

    /// \brief return the degeneracy
    int get_g() { return g; }

    /// \brief return the characteristic energy
    double get_omega() { return omega; }

    /// \brief evaluate the vibrational energy for a given quantum state
    double eval_E_vib( int iV );

    int get_Vmax() { return Vmax; }

private:
    int g;		// Degeneracy
    double omega;	// Characteristic energy

    int Vmax;
};

class Polyatom_electronic_level {
public:
    /// \brief Constructor
    Polyatom_electronic_level( std::vector<double> lev_data );
    
    /// \brief Deconstructor
    virtual ~Polyatom_electronic_level();
    
public:
    /// \brief Calculate maximum vibrational quantum numbers
    int calculate_V_max( int ivm );

    /// \brief Calculate the vibrational energy for a given set of quantum numbers
    double eval_E_vib( std::vector<int> iV);

    /// \brief Calculate rotational energy
    double eval_E_rot( int ivm, int iV, int iJ );

    /// \brief Calculate all the possible vibrational states
    void vib_loop( int im, std::vector<int> vib_state, std::vector< std::vector<int> > &results );

public:
    /* Fundamental level data */
    double N;
    double E;
    int g;
    
    /* Energy terms */
    double D;
    double r_e;
    
    /* Vibrational coefficients */
    std::vector<double> omega_e_vec;
    
    /* Vibrational modes */
    std::vector<Polyatom_vibrational_mode*> vmodes;

    /* Rotational characteristic frequencies */
    double A0;
    double B0;
    double C0;

    /* Rotational symmetry (?) factors */
    double sigma;
    double sigma_r;
};

#endif
