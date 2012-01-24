
#include <iostream>
#include <iomanip>

#include <string>
#include <sstream>

#include "flow_state.hh"

using namespace std;

Flow_state::
Flow_state() {}

Flow_state::
Flow_state(Gas_data &Q1, double u1, double Q_rad1)
 : u(u1), Q_rad(Q_rad1)
{
    Q = new Gas_data(Q1);
}

Flow_state::
Flow_state( Flow_state &fs )
 : u( fs.u ), Q_rad( fs.Q_rad )
{
    Q = new Gas_data(*fs.Q);
}

Flow_state::
~Flow_state() {}

void
Flow_state::
set_flow_state(Gas_data &Q1, double u1, double Q_rad1)
{
    Q = new Gas_data(Q1);
    u = u1;
    Q_rad = Q_rad1;
}

string
Flow_state::
str( bool with_Q_rad )
{
    ostringstream ost;
    ost << setprecision(12) << showpoint;
    for ( size_t itm=0; itm<Q->T.size(); ++ itm )
    	ost << setw(20) << Q->T[itm] << ' ';
    ost << setw(20) << Q->p << ' '
	<< setw(20) << Q->rho << ' '
	<< setw(20) << u << ' ';
    for( size_t isp=0; isp<Q->massf.size(); ++isp ) {
	ost << setw(20) << Q->massf[isp] << ' ';
    }
    if ( with_Q_rad ) {
    	ost << setw(20) << Q_rad << ' ';
    }

    return ost.str();

}
