// Author: Daniel F Potter
// Version: 13-Oct-2009
//          Initial coding.
// Version: 10-Dec-2009
//          Split into transport and diffusion flavours

#include <cmath>
#include <iostream>

#include "../../util/source/lua_service.hh"
#include "../../util/source/useful.h"
#include "../../nm/source/functor.hh"

#include "collision-integral.hh"
#include "chemical-species-library.hh"
#include "CI-functor.hh"
#include "physical_constants.hh"

using namespace std;

Collision_integral::Collision_integral( int iT, string type )
: iT_( iT ), type_( type )
{}

GuptaYos_CI_model::
GuptaYos_CI_model( int iT, int Z_i, int Z_j, lua_State * L )
: Collision_integral( iT, "GuptaYos_CI_model" ), Z_i_( Z_i ), Z_j_( Z_j )
{
    // Assume a table with the model paramaters is TOS.
    int ncurves = lua_objlen(L, -1);
    vector<Univariate_functor*> Pi_Omega_11;
    vector<Univariate_functor*> Pi_Omega_22;
    vector<Univariate_functor*> D;
    vector<double> breaks;
    double T_low, T_high;
    
    for ( int i = 0; i < ncurves; ++i ) {
    	lua_rawgeti(L, -1, i+1); // Lua indexes from 1 not 0
    	T_low = get_positive_number(L, -1, "T_low");
    	T_high = get_positive_number(L, -1, "T_high");
    	
    	Pi_Omega_11.push_back( new GuptaYos_CI_functor("Pi_Omega_11", T_low, T_high, L) );
    	Pi_Omega_22.push_back( new GuptaYos_CI_functor("Pi_Omega_22", T_low, T_high, L) );
    	D.push_back( new GuptaYos_CI_functor("D", T_low, T_high, L) );
    	breaks.push_back( T_low );
    	lua_pop(L,1);
    }
    breaks.push_back(T_high);
    
    Pi_Omega_11_ = new Segmented_functor( Pi_Omega_11, breaks );
    Pi_Omega_22_ = new Segmented_functor( Pi_Omega_22, breaks );
    D_ = new Segmented_functor( D, breaks );
    
    for ( int i = 0; i < ncurves; ++i ) {
	delete Pi_Omega_11[i];
	delete Pi_Omega_22[i];
	delete D[i];
    }
}

GuptaYos_CI_model::~GuptaYos_CI_model()
{
    delete Pi_Omega_11_;
    delete Pi_Omega_22_;
    delete D_;
}

double GuptaYos_CI_model::s_eval_Pi_Omega_11( Gas_data &Q )
{
    double Pi_Omega_11 = (*Pi_Omega_11_)(Q.T[iT_]) * 1.0e-20;	// convert Ang**2 -> m**2
    
    if ( Z_i_!=0 && Z_j_!=0 ) {
    	// if ( Q.p_e==0.0 ) {
    	//     cout << "GuptaYos_CI_model::s_eval_Pi_Omega_11()" << endl
    	//          << "Q.p_e = 0.0 -> this should not have passes the mole-fraction test" << endl;
    	// }
    	// We need to correct for the electron pressure
    	double tmpA = Q.T[iT_]/(1000.0*pow(Q.p_e/PC_P_atm,0.25));		// NOTE: using the electron temperature here
    	double f = 0.5*log(2.09e-2*pow(tmpA,4)+1.52*pow(tmpA,8.0/3.0));	// correction factor, Eq. 24b
    	Pi_Omega_11 *= f;
    }
    
    return Pi_Omega_11;	
}

double GuptaYos_CI_model::s_eval_Pi_Omega_22( Gas_data &Q )
{
    double Pi_Omega_22 = (*Pi_Omega_22_)(Q.T[iT_]) * 1.0e-20;	// convert Ang**2 -> m**2
    
    if ( Z_i_!=0 && Z_j_!=0 ) {
    	// if ( Q.p_e==0.0 ) {
    	//     cout << "GuptaYos_CI_model::s_eval_Pi_Omega_22()" << endl
    	//          << "Q.p_e = 0.0 -> this should not have passes the mole-fraction test" << endl;
    	// }
    	// We need to correct for the electron pressure
    	double tmpA = Q.T[iT_]/(1000.0*pow(Q.p_e/PC_P_atm,0.25));		// NOTE: using the electron temperature here
    	double f = 0.5*log(2.09e-2*pow(tmpA,4)+1.52*pow(tmpA,8.0/3.0));	// correction factor, Eq. 24b
    	Pi_Omega_22 *= f;
    }
    
    return Pi_Omega_22;	
}

double GuptaYos_CI_model::s_eval_D( Gas_data &Q )
{
    cout << "GuptaYos_CI_model::s_eval_D()" << endl
         << "WARNING: this function is obsolete and may not be returning the correct result." << endl;
    
    double D = (*D_)(Q.T[iT_]) / ( Q.p / PC_P_atm ) * 1.0e-4;	// convert cm**2-atm/sec - > m**2/sec
    
    if ( Z_i_!=0 && Z_j_!=0 ) {
    	if ( Q.p_e==0.0 ) {
    	    cout << "GuptaYos_CI_model::s_eval_D()" << endl
    	         << "Q.p_e = 0.0 -> this should not have passes the mole-fraction test" << endl;
    	}
    	// We need to correct for the electron pressure
    	double tmpA = Q.T[iT_]/(1000.0*pow(Q.p_e/PC_P_atm,0.25));		// NOTE: using the electron temperature here
    	double f = 2.0/log(2.09e-2*pow(tmpA,4)+1.52*pow(tmpA,8.0/3.0));		// correction factor, Eq. 42d
    	D *= f;
    }
    
    return D;
}

Stallcop_CI_model::
Stallcop_CI_model( int iT, int Z_i, int Z_j, lua_State * L )
: Collision_integral( iT, "Stallcop_CI_model" )
{
    // Free electron index
    ie_ = get_library_index_from_name( "e_minus" );
    
    // Create the Bivariate functors to calculate the collision-integrals
    // NOTE: no input file data required
    Pi_Omega_11_ = new Stallcop_SCP_CI_functor(1, 1, Z_i, Z_j);
    Pi_Omega_22_ = new Stallcop_SCP_CI_functor(2, 2, Z_i, Z_j);

}

Stallcop_CI_model::~Stallcop_CI_model()
{
    delete Pi_Omega_11_;
    delete Pi_Omega_22_;
}

double Stallcop_CI_model::s_eval_Pi_Omega_11( Gas_data &Q )
{
    // Electron number density
    double N_e = Q.massf[ie_] * Q.rho / PC_m_SI;
    
    double Pi_Omega_11 = (*Pi_Omega_11_)(Q.T[iT_],N_e) * 1.0e-20;	// convert Ang**2 -> m**2
    
    return Pi_Omega_11;	
}

double Stallcop_CI_model::s_eval_Pi_Omega_22( Gas_data &Q )
{
    // Electron number density
    double N_e = Q.massf[ie_] * Q.rho / PC_m_SI;
    
    double Pi_Omega_22 = (*Pi_Omega_22_)(Q.T[iT_],N_e) * 1.0e-20;	// convert Ang**2 -> m**2
    
    return Pi_Omega_22;	
}

double Stallcop_CI_model::s_eval_D( Gas_data &Q )
{
    cout << "Stallcop_CI_model::s_eval_D()" << endl
         << "WARNING: this function is obsolete and may not be returning the correct result." << endl;
    
    return 0.0;
}


Bruno_CI_model::
Bruno_CI_model( int iT, int Z_i, int Z_j, lua_State * L )
: Collision_integral( iT, "Bruno_CI_model" )
{
    // Create the Bivariate functors to calculate the collision-integrals
    if ( ( Z_i==0 && Z_j==0 ) || ( Z_i==1 && Z_j==0 ) || ( Z_i==0 && Z_j==1 ) ) {
    	// Heavy particle interaction
    	Pi_Omega_11_ = new Bruno_HPI_CI_functor(1, 1, Z_i, Z_j, L);
        Pi_Omega_22_ = new Bruno_HPI_CI_functor(2, 2, Z_i, Z_j, L);
    }
    else {
    	ostringstream oss;
    	oss << "Bruno_CI_model::Bruno_CI_model()" << endl
    	    << "Charge combination: Z_i = " << Z_i << ", Z_j = " << Z_j
    	    << " is not available." << endl;
    	input_error(oss);
    }
}

Bruno_CI_model::~Bruno_CI_model()
{
    delete Pi_Omega_11_;
    delete Pi_Omega_22_;
}

double Bruno_CI_model::s_eval_Pi_Omega_11( Gas_data &Q )
{
    return (*Pi_Omega_11_)(Q.T[iT_]) * 1.0e-20;				// convert Ang**2 -> m**2
}

double Bruno_CI_model::s_eval_Pi_Omega_22( Gas_data &Q )
{
    return (*Pi_Omega_22_)(Q.T[iT_]) * 1.0e-20;				// convert Ang**2 -> m**2
}

double Bruno_CI_model::s_eval_D( Gas_data &Q )
{
    cout << "Bruno_CI_model::s_eval_D()" << endl
         << "WARNING: this function is obsolete and may not be returning the correct result." << endl;
    
    return 0.0;
}


Neufeld_CI_model::
Neufeld_CI_model( int iT, Chemical_species * I, Chemical_species * J, lua_State * L )
: Collision_integral( iT, "Neufeld_CI_model" )
{
    // Extract the species charges
    int Z_i = I->get_Z();
    int Z_j = J->get_Z();
    
    // Create the univariate functors to calculate the collision-integrals
    if ( Z_i==0 && Z_j==0 ) {
    	// Heavy particle interaction
    	Pi_Omega_11_ = new Neufeld_CI_functor(1, 1, I, J, L);
        Pi_Omega_22_ = new Neufeld_CI_functor(2, 2, I, J, L);
    }
    else {
    	ostringstream oss;
    	oss << "Neufeld_CI_model::Neufeld_CI_model()" << endl
    	    << "Charge combination: Z_i = " << Z_i << ", Z_j = " << Z_j
    	    << " is not available." << endl;
    	input_error(oss);
    }
}

Neufeld_CI_model::~Neufeld_CI_model()
{
    delete Pi_Omega_11_;
    delete Pi_Omega_22_;
}

double Neufeld_CI_model::s_eval_Pi_Omega_11( Gas_data &Q )
{
    return (*Pi_Omega_11_)(Q.T[iT_]) * 1.0e-20;				// convert Ang**2 -> m**2
}

double Neufeld_CI_model::s_eval_Pi_Omega_22( Gas_data &Q )
{
    return (*Pi_Omega_22_)(Q.T[iT_]) * 1.0e-20;				// convert Ang**2 -> m**2
}

double Neufeld_CI_model::s_eval_D( Gas_data &Q )
{
    cout << "Neufeld_CI_model::s_eval_D()" << endl
         << "WARNING: this function is obsolete and may not be returning the correct result." << endl;
    
    return 0.0;
}

Collision_integral * new_CI_from_file( string fname )
{
    Collision_integral * CI_model = 0;

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);

    if( luaL_dofile(L, fname.c_str()) != 0 ) {
        ostringstream ost;
        ost << "new_CI_from_file():\n";
        ost << "Error in input file: " << fname << endl;
        input_error(ost);
    }

    lua_getglobal(L, "CI");
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "new_CI_from_file()\n";
        ost << "Expected a table entry for 'CI'" << endl;
        input_error(ost);
    }

    string i_name = get_string(L, -1, "i");
    string j_name = get_string(L, -1, "j");

    int Z_i = 0;
    int Z_j = 0;
    if ( i_name.find("_plus")!=string::npos ) Z_i = 1;
    else if ( i_name.find("_minus")!=string::npos ) Z_i = -1;
    if ( j_name.find("_plus")!=string::npos ) Z_j = 1;
    else if ( j_name.find("_minus")!=string::npos ) Z_j = -1;

    string CI_model_str = get_string(L,-1,"model");
    // Move to the parameters section
    lua_getfield(L, -1, "parameters" );
    if ( !lua_istable(L, -1) ) {
        ostringstream ost;
        ost << "new_CI_from_file()\n";
        ost << "Error setting parameters table for collision-integral: "
            << i_name << " - " << j_name << endl;
        input_error(ost);
    }

    // Create the CI model
    if ( CI_model_str=="none" ) {
        CI_model = new No_CI_model();
    }
    else if ( CI_model_str=="GuptaYos curve fits" ) {
        CI_model = new GuptaYos_CI_model( 0, Z_i, Z_j, L );
    }
    else if ( CI_model_str=="Stallcop curve fits" ) {
        CI_model = new Stallcop_CI_model( 0, Z_i, Z_j, L );
    }
    else if ( CI_model_str=="Bruno curve fits" ) {
        CI_model = new Bruno_CI_model( 0, Z_i, Z_j, L );
    }
    else if ( CI_model_str=="Neufeld curve fits" ) {
        ostringstream oss;
        oss << "new_CI_from_file()" << endl
            << "Collision pair: " << i_name << " - " << j_name << " requested the Neufeld curve fits model." << endl
            << "Not currently implementing this model." << endl;
        input_error( oss );
    }
    else {
        ostringstream oss;
        oss << "Binary_interaction::Binary_interaction()" << endl
            << "Collision integral model: " << CI_model_str << " is not recognised." << endl;
        input_error( oss );
    }
    lua_pop(L,1);       // pop 'parameters'

    lua_pop(L,1);       // pop 'CI'

    lua_close(L);

    return CI_model;
}

