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
    Polyatom_vibrational_mode( double omega, int Vmax, std::vector<int> Jmax );

    /// \brief Deconstructor
    virtual ~Polyatom_vibrational_mode() {};

public:
    double omega;

    int Vmax;
    std::vector<int> Jmax;
};

class Polyatom_electronic_level {
public:
    /// \brief Constructor
    Polyatom_electronic_level( std::vector<double> lev_data );
    
    /// \brief Deconstructor
    virtual ~Polyatom_electronic_level();
    
public:
    int calculate_V_max( int ivm );

    int calculate_J_max( int ivm, int iV );

    double eval_E_vib( int ivm, int iV );

    double eval_E_rot( int ivm, int iV, int iJ );

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

    /* Rotational coefficients */
    double A0;
    double B0;
    double C0;
    double sigma;
    double sigma_r;
};

#endif
