// Author: Rowan J. Gollan
// Date: 30-July-2008

#include <cmath>
#include <sstream>
#include <iostream>


#include "../../util/source/useful.h"
#include "../../util/source/lua_service.hh"
#include "hard-sphere-dcm.hh"
#include "physical_constants.hh"

using namespace std;

Hard_sphere_dcm::
Hard_sphere_dcm(lua_State *L)
    : R_mix_(0.0)
{
    lua_getglobal(L, "species");
    
    if ( !lua_istable(L, -1) ) {
	ostringstream ost;
	ost << "Hard_sphere_dcm::Hard_sphere_dcm():\n";
	ost << "Error in the declaration of species: a table is expected.\n";
	input_error(ost);
    }

    int nsp = lua_objlen(L, -1);
    d_.resize(nsp, 0.0);
    m_.resize(nsp, 0.0);
    R_.resize(nsp, 0.0);

    for ( int isp = 0; isp < nsp; ++isp ) {
	lua_rawgeti(L, -1, isp+1); // A Lua list is offset one from the C++ vector index
	const char* sp = luaL_checkstring(L, -1);
	lua_pop(L, 1);

	// Now bring the specific species table to TOS
	lua_getglobal(L, sp);
	if ( !lua_istable(L, -1) ) {
	    ostringstream ost;
	    ost << "Wilke_mixing_rule::Wilke_mixing_rule()\n";
	    ost << "Error locating information table for species: " << sp << endl;
	    input_error(ost);
	}

	d_[isp] = get_positive_value(L, -1, "d");
	
	// Get "R", store it, and then compute the mass of the particle.
	double M = get_positive_value(L, -1, "M");
	R_[isp] = PC_R_u/M;
	m_[isp] = M/PC_Avogadro; // to give kg

	lua_pop(L, 1); // Pops "sp" off stack
    }
    lua_pop(L, 1); // Pops "species" off stack
}

Hard_sphere_dcm::
~Hard_sphere_dcm() {}

int
Hard_sphere_dcm::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    
    // 0. update R
    R_mix_ = mass_average(Q.massf, R_);
    
    // 1. Proceed with calculation
    for( size_t i = 0; i < Q.D_AB.size(); ++i ) {
	for( size_t j = 0; j < Q.D_AB[i].size(); ++j ) {
	    if( i > j ) // already calculated
		Q.D_AB[i][j] = Q.D_AB[j][i];
	    Q.D_AB[i][j] = calculate_D_AB(Q, i, j);
	}
    }

    return SUCCESS;
}

double
Hard_sphere_dcm::
calculate_D_AB(Gas_data &Q, int i, int j)
{
    // Reference:
    // Bird, Stewart and Lightfoot (2001)
    // Transport Phenomena, 2nd edition
    // John Wiley & Sons, New York
    // Equation 17.3-10, p. 526

    const double small_n = 1.0e-15;

    double n_tot = Q.p / (PC_k_SI * Q.T[0]);
    double n_i = (Q.massf[i]*R_[i]/R_mix_)*n_tot; // first term gives mole fraction
    double n_j = (Q.massf[j]*R_[j]/R_mix_)*n_tot;
    double n = n_i + n_j;

    if( n <= small_n ) // At small_n, there's nothing to diffuse
	return 0.0;

    double d_ab = 0.5*(d_[i] + d_[j]);

    double tmpA = sqrt(PC_k_SI*Q.T[0]/M_PI);
    double tmpB = sqrt(0.5*((1.0/m_[i]) + (1.0/m_[j])));
    double tmpC = 1.0/(M_PI*d_ab*d_ab);
    double tmpD = 1.0/n;

    double D_AB = (2.0/3.0)*tmpA*tmpB*tmpC*tmpD;
    
    return D_AB;
}
