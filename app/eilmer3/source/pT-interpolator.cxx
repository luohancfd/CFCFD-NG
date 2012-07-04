// Author: Rowan J. Gollan
// Date: 15-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#include "../../../lib/util/source/useful.h"
#include "../../../lib/gas/models/gas-model.hh"
#include "cell.hh"
#include "kernel.hh"
#include "one_d_interp_scalar.hh"
#include "pT-interpolator.hh"

PT_interpolator::
PT_interpolator()
    : Thermo_interpolator() {}

PT_interpolator::
~PT_interpolator() {}

int
PT_interpolator::
s_one_d_interp(Gas_data &gL1, Gas_data &gL0,
	       Gas_data &gR0, Gas_data &gR1,
	       double cL1Length, double cL0Length,
	       double cR0Length, double cR1Length,
	       Gas_data &Lft, Gas_data &Rght)
{
    Gas_model *gmodel = get_gas_model_ptr();
    int apply_limiter_flag = get_apply_limiter_flag();
    int extrema_clipping_flag = get_extrema_clipping_flag();

    one_d_interp_scalar(gL1.p, gL0.p, gR0.p, gR1.p,
			cL1Length, cL0Length, cR0Length, cR1Length,
			Lft.p, Rght.p, apply_limiter_flag, extrema_clipping_flag);
    for ( int i = 0; i < gmodel->get_number_of_modes(); ++i ) {
	one_d_interp_scalar(gL1.T[i], gL0.T[i], gR0.T[i], gR1.T[i],
			    cL1Length, cL0Length, cR0Length, cR1Length,
			    Lft.T[i], Rght.T[i], apply_limiter_flag, extrema_clipping_flag);
    }

    int status = 0;

    if ( gmodel->eval_thermo_state_pT(Lft) != SUCCESS ) {
	printf("Rhoe_interpolator::s_one_d_interp(): Duff call to EOS.\n");
	printf("   Lft:\n");
	Lft.print_values();
	status = 1;
    }

    if ( gmodel->eval_thermo_state_pT(Rght) != SUCCESS ) {
	printf("Rhoe_interpolator::s_one_d_interp(): Duff call to EOS.\n");
	printf("   Rght:\n");
	Rght.print_values();
	status += 2;
    }

    if ( status == 0 )
	return SUCCESS;
    else
	return status;
}

