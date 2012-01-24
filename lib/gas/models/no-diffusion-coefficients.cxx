// Author: Rowan J. Gollan
// Date: 16-Jul-2008

#include "../../util/source/useful.h"
#include "no-diffusion-coefficients.hh"

using namespace std;

int
No_diffusion_coefficients::
s_eval_diffusion_coefficients(Gas_data &Q)
{
    for ( size_t i = 0; i < Q.D_AB.size(); ++i ) {
	for ( size_t j = 0; j < Q.D_AB[i].size(); ++j ) {
	    Q.D_AB[i][j] = 0.0;
	}
    }
    return SUCCESS;
}
