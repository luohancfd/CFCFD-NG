// Author: Daniel F Potter
// Version: 13-Apr-2010
//          Initial coding.

#include <cmath>
#include <iostream>

#include "../../util/source/lua_service.hh"

#include "physical_constants.hh"
#include "chemical-species-library.hh"
#include "binary-interaction.hh"

using namespace std;

Binary_interaction::
Binary_interaction( int isp, int jsp, lua_State *L )
: isp_( isp ), jsp_( jsp )
{
    // use chemical species library to get molecular weights and T index
    if ( ! chemical_species_library_initialised() ) {
    	ostringstream ost;
	ost << "Binary_interaction::Binary_interaction()\n";
	ost << "The chemical species library must be initialised.\n";
	input_error(ost);
    }
    Chemical_species * I = get_library_species_pointer(isp);
    Chemical_species * J = get_library_species_pointer(jsp);
    m_i_ = I->get_M() / PC_Avogadro;
    m_j_ = J->get_M() / PC_Avogadro;
    if ( I->get_name()=="e_minus" ) iT_ = I->get_iT_trans();
    else if ( J->get_name()=="e_minus" ) iT_ = J->get_iT_trans();
    else iT_ = I->get_iT_trans();	// should be the same as for species 'J'
    
    // Get the collision integral model type
    string CI_model_str = get_string(L,-1,"model");
    // Move to the parameters section
    lua_getfield(L, -1, "parameters" );
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Binary_interaction::Binary_interaction()\n";
	ost << "Error setting parameters table for collision-integral: " 
	    << I->get_name() << " - " << J->get_name() << endl;
	input_error(ost);
    }

    // Create the CI model
    if ( CI_model_str=="none" ) {
    	CI_model_ = new No_CI_model();
    }
    else if ( CI_model_str=="GuptaYos curve fits" ) {
    	CI_model_ = new GuptaYos_CI_model( iT_, I->get_Z(), J->get_Z(), L );
    }
    else if ( CI_model_str=="Stallcop curve fits" ) {
    	CI_model_ = new Stallcop_CI_model( iT_, I->get_Z(), J->get_Z(), L );
    }
    else if ( CI_model_str=="Bruno curve fits" ) {
    	CI_model_ = new Bruno_CI_model( iT_, I->get_Z(), J->get_Z(), L );
    }
    else if ( CI_model_str=="Neufeld curve fits" ) {
    	CI_model_ = new Neufeld_CI_model( iT_, I, J, L );
    }
    else {
	ostringstream oss;
	oss << "Binary_interaction::Binary_interaction()" << endl
	    << "Collision integral model: " << CI_model_str << " is not recognised." << endl;
	input_error( oss );
    }
    lua_pop(L,1);	// pop 'parameters'
}

Binary_interaction::
~Binary_interaction()
{
    delete CI_model_;
}

double
Binary_interaction::
s_eval_Delta_1( Gas_data &Q )
{
    // NOTE: all parameters now in SI units
    double Pi_Omega_11 = CI_model_->eval_Pi_Omega_11( Q );
    // if ( !finite(Pi_Omega_11) ) {
    	// cout << "type = " << CI_model_->get_type() << endl;
    // }
    return (8.0/3.0)*(pow((2.0*m_i_*m_j_)/(M_PI*PC_k_SI*Q.T[iT_]*(m_i_+m_j_)),0.5))*Pi_Omega_11;
}

double
Binary_interaction::
s_eval_Delta_2( Gas_data &Q )
{
    // NOTE: all parameters now in SI units
    double Pi_Omega_22 = CI_model_->eval_Pi_Omega_22( Q );
    // if ( !finite(Pi_Omega_22) ) {
    	// cout << "type = " << CI_model_->get_type() << endl;
    // }
    return (16.0/5.0)*(pow((2.0*m_i_*m_j_)/(M_PI*PC_k_SI*Q.T[iT_]*(m_i_+m_j_)),0.5))*Pi_Omega_22;
}

double
Binary_interaction::
s_eval_D( Gas_data &Q )
{
    // NOTE: Evaluating from definition rather than curve-fit
    double Delta_1 = s_eval_Delta_1(Q);
    double D_ij = PC_k_SI * Q.T[iT_] / ( Q.p * Delta_1 );
    // cout << "D_ij (calc) = " << D_ij << endl;
    // cout << "D_ij (CI_model) = " << CI_model_->eval_D( Q ) << endl;
    
    return D_ij;
    // return CI_model_->eval_D( Q );
}

