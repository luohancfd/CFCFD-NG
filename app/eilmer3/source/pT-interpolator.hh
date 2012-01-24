// Author: Rowan J. Gollan
// Date: 15-May-2009
// Place: NASA Langley, Hampton, Virginia, USA
//

#ifndef PT_INTERPOLATOR_HH
#define PT_INTERPOLATOR_HH

#include "thermo-interpolator.hh"

class PT_interpolator : public Thermo_interpolator {
public:
    PT_interpolator();
    ~PT_interpolator();

private:
    int s_one_d_interp(Gas_data &gL1, Gas_data &gL0,
		       Gas_data &gR0, Gas_data &gR1,
		       double cL1Length, double cL0Length,
		       double cR0Length, double cR1Length,
		       Gas_data &Lft, Gas_data &Rght);
};

#endif
