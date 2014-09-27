// Author: Daniel. F Potter
// Date: 24-Sep-2009

#include "../../util/source/useful.h"
#include "no-transport-coefficients.hh"

using namespace std;

int
No_transport_coefficients::
s_eval_transport_coefficients(Gas_data &Q, Gas_model *gmodel)
{
    // 1. Viscosity coefficient
    Q.mu = 0.0;
    
    // 2. Thermal conductivities
    for ( size_t itm = 0; itm < Q.k.size(); ++itm ) {
    	Q.k[itm] = 0.0;
    }
	
    return SUCCESS;
}

