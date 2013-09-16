/** \file stark_width_models.cxx
 *  \ingroup radiation
 *
 *  \brief Classes for describing Stark broadening of lines
 *
 *  \author Daniel F. Potter
 *  \version 16-Sept-13: initial implementation
 *            
 **/

#include <cstdlib>
#include <string>
#include <math.h>

#include "stark_width_models.hh"
#include "radiation_constants.hh"

using namespace std;

// Base class

StarkWidthModel::StarkWidthModel( string name, double n )
 : name( name ), n( n ) {}

StarkWidthModel::~StarkWidthModel()
{}

// Approximate model using a curve fit for gamma_s0

ApproxStarkWidth::ApproxStarkWidth( double n, double constA, double constB,
                                              double      I, double    E_u )
 : StarkWidthModel("ApproxStarkWidth",n), constA( constA ), constB( constB ),
   I_icm( I/(RC_c*RC_h_SI) ), E_u_icm( E_u/(RC_c*RC_h_SI) )
{}

ApproxStarkWidth::~ApproxStarkWidth()
{}

double ApproxStarkWidth::eval( double T_e, double N_e )
{
    // NOTE: applying fabs() around delta_E as E_u can be > E_ionise
    double tmp =  0.5 * constA / pow( fabs(I_icm - E_u_icm), constB );
    double gamma_S0 = tmp * RC_c;
    double gamma_S = gamma_S0 * pow( (T_e / 1.0e4), n ) * ( N_e / 1.0e16 );

    return gamma_S;
}

// 'Exact' model using reference width and exponential factor obtained from
// curve fits to the tabulated data of Griem

GriemStarkWidth::GriemStarkWidth( double n, double gamma_S0 )
 : StarkWidthModel("GriemStarkWidth",n), gamma_S0( gamma_S0 )
{}

GriemStarkWidth::~GriemStarkWidth()
{}

double GriemStarkWidth::eval( double T_e, double N_e )
{
    // NOTE: applying fabs() around delta_E as E_u can be > E_ionise
    double gamma_S = gamma_S0 * pow( (T_e / 1.0e4), n ) * ( N_e / 1.0e16 );

    return gamma_S;
}
