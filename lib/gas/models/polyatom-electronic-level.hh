// Author: Daniel F. Potter
// Version: 27-July-2012
//          Ported from lib/gas/models/diatom-electronic-level.hh

#ifndef POLYATOM_ELEC_LEV_HH
#define POLYATOM_ELEC_LEV_HH

#include <vector>

#include "gas_data.hh"

#define FULL_ROVIBRATIONAL_SUMMATION 0

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

    /// \brief evaluate the vibrational partition function
    double eval_Q_from_T( double T );

    /// \brief return the maximum vibrational quantum number
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
    virtual double eval_E_rot( int iJ, int iK ) = 0;

    /// \brief Calculate all the possible vibrational states
    void vib_loop( int im, std::vector<int> &vib_state, std::vector< std::vector<int> > &results );

    /// \brief Calculate the equilibrium partition function for this level
    virtual double eval_Q_from_T( double T) = 0;

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

class Spherical_top_polyatom_electronic_level : public Polyatom_electronic_level {
public:
    /// \brief Constructor
    Spherical_top_polyatom_electronic_level( std::vector<double> lev_data );

    /// \brief Deconstructor
    ~Spherical_top_polyatom_electronic_level();

public:
    /// \brief Calculate rotational energy
    double eval_E_rot( int iJ, int iK );

    /// \brief Calculate the equilibrium partition function for this level
    double eval_Q_from_T( double T);
};

class Symmetrical_top_polyatom_electronic_level : public Polyatom_electronic_level {
public:
    /// \brief Constructor
    Symmetrical_top_polyatom_electronic_level( std::vector<double> lev_data );

    /// \brief Deconstructor
    ~Symmetrical_top_polyatom_electronic_level();

public:
    /// \brief Calculate rotational energy
    double eval_E_rot( int iJ, int iK );

    /// \brief Calculate the equilibrium partition function for this level
    double eval_Q_from_T( double T);
};

class Asymmetric_top_polyatom_electronic_level : public Polyatom_electronic_level {
public:
    /// \brief Constructor
    Asymmetric_top_polyatom_electronic_level( std::vector<double> lev_data );

    /// \brief Deconstructor
    ~Asymmetric_top_polyatom_electronic_level();

public:
    /// \brief Calculate rotational energy
    double eval_E_rot( int iJ, int iK );

    /// \brief Calculate the equilibrium partition function for this level
    double eval_Q_from_T( double T);
};

#endif
